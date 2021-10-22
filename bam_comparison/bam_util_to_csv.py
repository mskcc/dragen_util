#!/usr/bin/env python

import sys
import os
import pysam

USAGE="python bam_util_to_csv.py control.bam experimental.bam"

# Save all errors until end and then log
ERRORS = []
TOTAL_MISSING_READS = 0

# Constants to track value differences between the two BAMs
CONTROL_BAM = 'CONTROL_BAM'
TARGET_BAM = 'TARGET_BAM'

# Confidence levels for pairing - the more values that are the same, the higher confidence. The interesting pairs are
# medium/low-confidence matches (CONFIDENCE_MEDIUM/LOW) b/c for these there is a difference we would care about
CONFIDENCE_HIGH = 'high'            # High-Confidence Match:    [ read_id, flag, mapq ]
CONFIDENCE_MEDIUM = 'medium'        # Medium-Confidence Match:  [ read_id, flag ]
CONFIDENCE_LOW = 'low'              # Low-Confidence Match:     [ read_id ]

def combine_flag_vals(dic, flag_list, val_type):
    """Extends all the values across all flags of an Entry:flag_to_val_dic field

    :param Entry:flag_to_val_dic dic:   flag-to-entries dictionary
    :param int[] flag_list:             list of flags, e.g. [ 63, 1023,... ]
    :param string, val_type:            CONTROL_BAM/TARGET_BAM
    :return: int[ int[] ]               list of [flag, value] pairs
    """
    flag_read_info_entries = []
    for f in flag_list:
        read_info_entries = dic[f][val_type]
        flag_read_info_entries.extend(read_info_entries)

    return flag_read_info_entries

def record_unpaired_read(read_info_list, read_info_type, bam_type):
    """ Record unpaired reads in the input read_info list
            - Adds to global TOTAL_MISSING_READS
            - Adds to global ERRORS list
        Output csv list:
            e.g. [ "UNPAIRED_FW,CONTROL_BAM,read_id,reference_name,flag,mapq",... ]
    """
    global TOTAL_MISSING_READS
    global ERRORS

    for read_info_entry in read_info_list:
        TOTAL_MISSING_READS += 1
        ERRORS.append("{},{},{},{},{},{}".format(
            read_info_type,
            bam_type,
            read_info_entry.get_read(),
            read_info_entry.get_refr(),
            str(read_info_entry.get_flag()),
            str(read_info_entry.get_mapq())
        ))

class Entry:
    """ Tracks all the records for an input read name (SAM QNAME)

    :field read:                        QNAME of SAM record
    :field flag_to_val_dic:             flag -> { CONTROL/TARGET -> [ MAPQ_values ] }
    """
    read = None
    flag_to_val_dic = {}

    def __init__(self, read):
        self.read = read
        self.flag_to_val_dic = {}

    def add_flag_paths(self, flag):
        """ Safely adds path for a flag (populates if the path in the dictionary doesn't exist)
        """
        if flag not in self.flag_to_val_dic:
             self.flag_to_val_dic[flag] = {}
        val_dic = self.flag_to_val_dic[flag]
        if CONTROL_BAM not in val_dic:
            val_dic[CONTROL_BAM] = []
        if TARGET_BAM not in val_dic:
            val_dic[TARGET_BAM] = []

    def add_v1(self, read_info):
        """ Adds the CONTROL_BAM flag-MAPQ_value pair
        """
        flag = read_info.get_flag()
        self.add_flag_paths(flag)
        self.flag_to_val_dic[flag][CONTROL_BAM].append(read_info)

    def add_v2(self, read_info):
        """ Adds the TARGET flag-MAPQ_value pair
        """
        flag = read_info.get_flag()
        self.add_flag_paths(flag)
        self.flag_to_val_dic[flag][TARGET_BAM].append(read_info)

    def return_flag_val_pairs(self):
        # Note - this permanently modifies the flag_to_val_dic field
        pairs = []

        paired_up_flags = []
        # Add pairs

        for flag, val_dic in self.flag_to_val_dic.items():
            control_reads = val_dic[CONTROL_BAM]
            target_reads = val_dic[TARGET_BAM]

            # Dic, { mapq -> Read_Info[] }, for (read_id, flag) pair that maps CONTROL MAPQ scores to their read info
            # Later, the TARGET read_info instances will be mapped based on matching mapq scores
            read_info_pairing_dic = {}
            while len(control_reads) > 0:
                control_read_info = control_reads.pop()
                control_mapq = control_read_info.get_mapq()
                if control_mapq in read_info_pairing_dic:
                    read_info_pairing_dic[control_mapq].append(control_read_info)
                else:
                    read_info_pairing_dic[control_mapq] = [ control_read_info ]

            unpaired_target_reads = []
            while len(target_reads) > 0:
                target_read_info = target_reads.pop()
                target_mapq = target_read_info.get_mapq()
                if target_mapq in read_info_pairing_dic:
                    # High-Confidence Match: (read_id, flag, mapq)
                    matched_control_read_info = read_info_pairing_dic[target_mapq].pop()
                    pairs.append([ CONFIDENCE_HIGH, matched_control_read_info, target_read_info ])        # Same value

                    # Remove if entry for mapq if all control_bams of that (read_id, flag) pair have been matched
                    if len(read_info_pairing_dic[target_mapq]) == 0:
                        del read_info_pairing_dic[target_mapq]
                else:
                    unpaired_target_reads.append(target_read_info)

            # Extract all unpaired control Read_Info's to a list
            unpaired_control_reads = []
            for mapq, read_info_list in read_info_pairing_dic.items():
                 unpaired_control_reads.extend(read_info_list)

            # Pair up any remaining reads between target & control
            # TODO - We would probably want to pair reads that have the closest scores b/c the MAPQ differences could
            # vary significantly if there are multiple control/target options w/ different MAPQ scores
            # Right now, we approximate this w/ sorting...
            sorted_unpaired_control_reads = sorted(unpaired_control_reads, key=lambda read_info: read_info.get_mapq())
            sorted_unpaired_target_reads = sorted(unpaired_target_reads, key=lambda read_info: read_info.get_mapq())
            while len(sorted_unpaired_control_reads) > 0 and len(sorted_unpaired_target_reads) > 0:
                # Medium-Confidence Match: [ (read_id, flag) ]
                pairs.append( [ CONFIDENCE_MEDIUM, sorted_unpaired_control_reads.pop(), sorted_unpaired_target_reads.pop()] )

            if len(sorted_unpaired_control_reads) == 0 and len(sorted_unpaired_target_reads) == 0:
                # All values for flag are paired - keep track of these flags so they can be deleted later
                paired_up_flags.append(flag)
            else:
                # Re-assign unpaired values for flag
                self.flag_to_val_dic[flag][CONTROL_BAM] = sorted_unpaired_control_reads
                self.flag_to_val_dic[flag][TARGET_BAM] = sorted_unpaired_target_reads

        # Remove flags for a read that have been successfully paired
        for flag in paired_up_flags:
            del self.flag_to_val_dic[flag]

        # Try to pair flags that remain w/ only control/target values

        # Separate flags by reverse-complement/forward & target/control
        remaining_flags = self.flag_to_val_dic.keys()
        rc_control_flags = [ f for f in remaining_flags if is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][CONTROL_BAM]) > 0 ]
        rc_target_flags = [ f for f in remaining_flags if is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][TARGET_BAM]) > 0 ]
        fw_control_flags = [ f for f in remaining_flags if not is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][CONTROL_BAM]) > 0 ]
        fw_target_flags = [ f for f in remaining_flags if not is_flag_reverse_complemented(f) and len(self.flag_to_val_dic[f][TARGET_BAM]) > 0 ]

        # Combine read_info entries across flags
        rc_control_vals = combine_flag_vals(self.flag_to_val_dic, rc_control_flags, CONTROL_BAM)
        rc_target_vals = combine_flag_vals(self.flag_to_val_dic, rc_target_flags, TARGET_BAM)
        fw_control_vals = combine_flag_vals(self.flag_to_val_dic, fw_control_flags, CONTROL_BAM)
        fw_target_vals = combine_flag_vals(self.flag_to_val_dic, fw_target_flags, TARGET_BAM)

        # Pair read_info entries
        while len(rc_control_vals) > 0 and len(rc_target_vals) > 0:
            rc_control_read_info = rc_control_vals.pop()
            rc_target_read_info = rc_target_vals.pop()
            # Low-Confidence Match: [ read_id ]
            pairs.append([ CONFIDENCE_LOW, rc_control_read_info, rc_target_read_info ])
        while len(fw_control_vals) > 0 and len(fw_target_vals) > 0:
            fw_control_read_info = fw_control_vals.pop()
            fw_target_read_info = fw_target_vals.pop()
            # Low-Confidence Match: [ read_id ]
            pairs.append([ CONFIDENCE_LOW, fw_control_read_info, fw_target_read_info ])

        record_unpaired_read(rc_control_vals, "UNPAIRED_RC", CONTROL_BAM)
        record_unpaired_read(rc_target_vals, "UNPAIRED_RC", TARGET_BAM)
        record_unpaired_read(fw_control_vals, "UNPAIRED_FW", CONTROL_BAM)
        record_unpaired_read(fw_target_vals, "UNPAIRED_FW", TARGET_BAM)

        return pairs

class Read_Info():
    """ Simple record of the read information, excludes qname/read_id b/c it's accounted for elsewhere

    :field score
    :field flag
    :field refr     scaffold read aligned to
    """
    score = None
    flag = None
    refr = None
    read = None
    def __init__(self, score, flag, refr, read):
        self.score = score
        self.flag = flag
        self.refr = refr
        self.read = read

    def get_mapq(self):
        return self.score

    def get_flag(self):
        return self.flag

    def get_refr(self):
        return self.refr

    def get_read(self):
        return self.read

    def to_string(self):
        return "{},{},{}" % (self.read, int(self.flag), self.refr)

def is_flag_reverse_complemented(flag):
    # 5th bit indicates reverse complement
    fifth_bit = get_kth_bit(flag, 5)
    return fifth_bit == 1

def get_kth_bit(n, k):
    """ Returns kth bit of an int

    REF: https://www.geeksforgeeks.org/find-value-k-th-bit-binary-representation/
    """
    return (n & (1 << (k - 1))) >> (k - 1)

def fail(err_msg = None):
    """ Exits Application and logs error message """
    print("Usage: %s" % USAGE)
    if(err_msg):
        print("ERROR: " + err_msg)
    sys.exit(1)

def write_file(file_name, entry_list):
    """ Writes @contents to @file_name
    :param file_name:
    :param contents:
    :return:
    """
    merge_commands_file = open(file_name, "a")
    merge_commands_file.truncate(0)  # Delete any old data
    merge_commands_file.write("\n".join(entry_list))
    merge_commands_file.close()

def extract_bam_entries(bam_file, entry_map, type):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    aln_segments = samfile.fetch()
    for aln_seg in aln_segments:
        read_id = aln_seg.query_name
        flag = aln_seg.flag
        score = aln_seg.mapping_quality
        refr = aln_seg.reference_name

        if read_id not in entry_map:
            entry_map[read_id] = Entry(read_id)
        entry = entry_map[read_id]

        ri_entry = Read_Info(score, flag, refr, read_id)
        if type == CONTROL_BAM:
            entry.add_v1(ri_entry)
        elif type == TARGET_BAM:
            entry.add_v2(ri_entry)

    return entry_map

def parse_entries(b1, b2):
    entry_map = {}
    entry_map = extract_bam_entries(b1, entry_map, CONTROL_BAM)
    entry_map = extract_bam_entries(b2, entry_map, TARGET_BAM)

    entries = []
    for read_id, entry in entry_map.items():
        control_target_read_info_pairs = entry.return_flag_val_pairs()
        for ct_pair in control_target_read_info_pairs:
            confidence_level = ct_pair[0]
            control_read_info = ct_pair[1]
            target_read_info = ct_pair[2]

            control_flag = control_read_info.get_flag()
            target_flag = target_read_info.get_flag()

            control_mapq = control_read_info.get_mapq()
            target_mapq = target_read_info.get_mapq()

            control_refr = control_read_info.get_refr()
            target_refr = target_read_info.get_refr()

            line = "{},{},{},{},{},{},{},{},{},{}".format(
                read_id,
                confidence_level,
                control_refr,
                target_refr,
                str(control_flag),
                str(target_flag),
                str(control_mapq),
                str(target_mapq),
                control_mapq - target_mapq,
                target_mapq - control_mapq
            )
            entries.append(line)

    return entries

def summarize(paired_reads):
    high_confidence = [ line for line in paired_reads if ",{},".format(CONFIDENCE_HIGH) in line ]
    med_confidence = [ line for line in paired_reads if ",{},".format(CONFIDENCE_MEDIUM) in line ]
    low_confidence = [ line for line in paired_reads if ",{},".format(CONFIDENCE_LOW) in line ]
    print("Number Pair(s): %d" % len(paired_reads))
    print("\t High-Confidence Pair(s): %d" % len(high_confidence))
    print("\t Medium-Confidence Pair(s): %d" % len(med_confidence))
    print("\t Low-Confidence Pair(s): %d" % len(low_confidence))
    print("Unpaired Read(s): %d" % TOTAL_MISSING_READS)

def main():
    if len(sys.argv) != 3:
        fail("Please specify two bam files")

    inputs = sys.argv[1:]
    for input in inputs:
        if not os.path.isfile(input):
            fail("%s is not a valid file" % input)

    b1 = inputs[0]
    b2 = inputs[1]

    basename = "{}____{}".format(b1.split("/")[-1].split(".")[0], b2.split("/")[-1].split(".")[0])
    output_file = "{}___bam_differences.csv".format(basename)
    missing_file = "{}___missing.csv".format(basename)

    print("%s=%s\n%s=%s\nOUTPUT=%s" % (CONTROL_BAM, b1, TARGET_BAM, b2, output_file))

    paired_reads = parse_entries(b1, b2)

    print("Writing paired CSV: %s" % output_file)
    paired_content = [ "qname,confidence,control_refr,target_refr,control_flag,target_flag,control_mapq,target_mapq,ct_diff,tc_diff"]
    paired_content.extend(paired_reads)
    write_file(output_file, paired_content)

    print("Writing unpaired CSV: %s" % missing_file)
    unpaired_content = [ "read_info_type,bam,qname,refr,flag,mapq" ]
    unpaired_content.extend(ERRORS)
    write_file(missing_file, unpaired_content)

    summarize(paired_reads)

if __name__ == '__main__':
    main()
