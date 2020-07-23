# Test applymask.py
# Run all tests: python3 test_applymask.py
# Run one test:  python3 test_applymask.py TestApplyMask.test_load_mask_range_aux
# Run code coverage: coverge run test_applymask.py 
# View code coverage report: coverage report -m
# Generate code coverage html report: coverage html

import os
import unittest
import applymask

def assert_fasta(self, fasta_name, sequence):
    output_masked_fasta = fasta_name
    header, seq, chunklen = applymask.load_fasta(output_masked_fasta)
    self.assertEqual(header, self.header)
    self.assertEqual(chunklen, self.chunklen)
    self.assertEqual(seq, sequence)           

def assert_fasta_zip(self, zip_fasta_name, sequence):
    output_masked_fasta_gz = zip_fasta_name
    header, seq, chunklen = applymask.load_fasta_gzip(output_masked_fasta_gz)
    self.assertEqual(header, self.header)
    self.assertEqual(chunklen, self.chunklen)
    self.assertEqual(seq, sequence)        

class TestApplyMask(unittest.TestCase):
    def setUp(self):
        self.header = ">NC_000962_3"
        self.chunklen = 60
        self.expected_mask= "".join([
                            "100000000000000000000000000000000000000000000000000000000001",
                            "100000000000000000000000000000000000000000000000000000000001", 
                            "100000000000000000000000000000000000000000000000000000000001", 
                            "100000000000000000000000000000000000000000000000000000000001", 
                            "100000000000000000000000000000000000000000000000000000000001", 
                            "000000000000000000000000000000000000000000000000000000000000", 
                            "000000000000000000000000000000000000000000000000000000000000", 
                            "000000000000000000000000000000000000000000000000000000000000", 
                            "000000000000000000000000000000000000000000000000000000000000", 
                            "000000000000000000000000000000000000000000000000000000000000" ])
        self.sequence = "".join([
                            "TTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTT",
                            "AACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTG", 
                            "ACCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTT", 
                            "GCTCTGTTATCCGTGCCGAGCAGCTTTGTCCAAAACGAAATCGAGCGCCATCTGCGGGCC", 
                            "CCGATTACCGACGCTCTCAGCCGCCGACTCGGACATCAGATCCAACTCGGGGTCCGCATC", 
                            "GCTCCGCCGGCGACCGACGAAGCCGACGACACTACCGTGCCGCCTTCCGAAAATCCTGCT", 
                            "ACCACATCGCCAGACACCACAACCGACAACGACGAGATTGATGACAGCGCTGCGGCACGG", 
                            "GGCGATAACCAGCACAGTTGGCCAAGTTACTTCACCGAGCGCCCGCACAATACCGATTCC", 
                            "GCTACCGCTGGCGTAACCAGCCTTAACCGTCGCTACACCTTTGATACGTTCGTTATCGGC", 
                            "GCCTCCAACCGGTTCGCGCACGCCGCCGCCTTGGCGATCGCAGAAGCACCCGCCCGCGCT" ])
        self.maskedsequence = "".join([
                            "NTGACCGATGACCCCGGTTCAGGCTTCACCACAGTGTGGAACGCGGTCGTCTCCGAACTN",
                            "NACGGCGACCCTAAGGTTGACGACGGACCCAGCAGTGATGCTAATCTCAGCGCTCCGCTN", 
                            "NCCCCTCAGCAAAGGGCTTGGCTCAATCTCGTCCAGCCATTGACCATCGTCGAGGGGTTN", 
                            "NCTCTGTTATCCGTGCCGAGCAGCTTTGTCCAAAACGAAATCGAGCGCCATCTGCGGGCN", 
                            "NCGATTACCGACGCTCTCAGCCGCCGACTCGGACATCAGATCCAACTCGGGGTCCGCATN", 
                            "GCTCCGCCGGCGACCGACGAAGCCGACGACACTACCGTGCCGCCTTCCGAAAATCCTGCT", 
                            "ACCACATCGCCAGACACCACAACCGACAACGACGAGATTGATGACAGCGCTGCGGCACGG", 
                            "GGCGATAACCAGCACAGTTGGCCAAGTTACTTCACCGAGCGCCCGCACAATACCGATTCC", 
                            "GCTACCGCTGGCGTAACCAGCCTTAACCGTCGCTACACCTTTGATACGTTCGTTATCGGC", 
                            "GCCTCCAACCGGTTCGCGCACGCCGCCGCCTTGGCGATCGCAGAAGCACCCGCCCGCGCT" ])
    @classmethod
    def tearDownClass(cls):
        os.remove("test.masked.fasta")
        os.remove("test.fasta.masked.gz")

    def test_list_to_range_str_1(self):
        data = [1,2,3,4,5,7]
        expected = [(1,5),(7,7)]

        result = applymask.lst_to_range_str(data)
        self.assertEqual(expected, result)
    
    def test_list_to_range_str_2(self):
        data = [0, 59, 60, 119, 120, 179, 180, 239, 240, 299]
        expected = [(0,0),(59,60),(119,120), (179,180), (239,240), (299,299)]

        result = applymask.lst_to_range_str(data)
        self.assertEqual(expected, result)

    def test_list_to_range_str_3(self):
        data = [10, 59, 60, 119, 120, 179, 180, 239, 240, 299]
        expected = [(10,10),(59,60),(119,120), (179,180), (239,240), (299,299)]

        result = applymask.lst_to_range_str(data)
        self.assertEqual(expected, result)

    def test_load_mask_range_aux(self):
        data = ["0\t0","5\t7", "9\t9"]
        expected = list("1000011101")

        result = applymask.load_mask_range_aux(data, 10)
        self.assertEqual(result, expected)

    def test_load_mask_range_aux_exception(self):
        data = ["0\ta","5\t7", "9\t9"]
        self.assertRaises(Exception, applymask.load_mask_range_aux(data,10))

    def test_load_mask_position_aux(self):
        data = ["0","5", "6", "7", "9"]
        expected = list("1000011101")

        result = applymask.load_mask_position_aux(data, 10)
        self.assertEqual(result, expected)

    def test_load_mask_position_aux_exception(self):
        data = ["a","5", "6", "7", "9"]
        self.assertRaises(Exception, applymask.load_mask_position_aux(data,10))
    
    def test_get_mask_ranges(self):
        data = list("1000011101")
        expected = "\n".join(["0\t0","5\t7", "9\t9"])
        
        result = applymask.get_mask_ranges(data)
        self.assertEqual(result, expected)

    def test_apply_mask(self):
        mask = list("1000011101")
        sequence = list("ACGTACGTAC")

        expected = "NCGTANNNAN"
        result = applymask.apply_mask(mask, sequence)
        self.assertEqual(result, expected)

    def test_load_mask_fasta(self):
        inputfile = "data/mask_fasta.fasta"
        expected = self.expected_mask

        result = applymask.load_mask_fasta(inputfile)
        self.assertEqual(result, expected)
    
    def test_load_mask_position(self):
        inputfile = "data/mask_position.txt"
        expected = list(self.expected_mask)

        result = applymask.load_mask_positon(inputfile,600)
        self.assertEqual(result, expected)

    def test_load_mask_range(self):
        inputfile = "data/mask_range.tsv"
        expected = list(self.expected_mask)

        result = applymask.load_mask_range(inputfile,600)
        self.assertEqual(result, expected)

    def test_load_fasta(self):
        assert_fasta(self,"data/test.fasta", self.sequence)
    
    def test_load_fasta_gzip(self):
        inputfile = "data/test.fasta.gz"
        expected_sequence = self.sequence
        assert_fasta_zip(self, inputfile, expected_sequence)


    def test_string_insert_newlines(self):
        input_str = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
        expected_str = "AAAAAAAAAA\nCCCCCCCCCC\nGGGGGGGGGG\nTTTTTTTTTT"

        result = applymask.string_insert_newlines(input_str,10)
        self.assertEqual(expected_str, result)

    def test_save_fasta(self):
        filepath= "data/saved_fasta.fasta"
        header = ">NC_000962_3"
        sequence = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
        chunklen = 10
        applymask.save_fasta(filepath,header,sequence,chunklen)
        header_out, seq_out, chunklen_out = applymask.load_fasta(filepath)
        self.assertEqual(header, header_out)
        self.assertEqual(chunklen, chunklen_out)
        self.assertEqual(sequence,seq_out)
        os.remove(filepath)
        print(f"file {filepath} removed.")

    def test_save_fasta_gzip(self):
        filepath= "data/saved_fasta.fasta.gz"
        header = ">NC_000962_3"
        sequence = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT"
        chunklen = 10
        applymask.save_fasta_gzip(filepath,header,sequence,chunklen)
        header_out, seq_out, chunklen_out = applymask.load_fasta_gzip(filepath)
        self.assertEqual(header, header_out)
        self.assertEqual(chunklen, chunklen_out)
        self.assertEqual(sequence,seq_out)
        os.remove(filepath)
        print(f"file {filepath} removed.")
    
    def test_main_fasta_nozip(self):
        fasta_filepath = "data/test.fasta"
        mask_filepath = "data/mask_fasta.fasta"
        mask_format = "fasta"
        use_gzip = "false"
        print_mask_ranges = "true"
        output_fasta = "test.masked.fasta"
        applymask.main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges)
        assert_fasta(self, output_fasta, self.maskedsequence)

    def test_main_position_nozip(self):
        fasta_filepath = "data/test.fasta"
        mask_filepath = "data/mask_position.txt"
        mask_format = "position"
        use_gzip = "false"
        print_mask_ranges = "true"
        output_fasta = "test.masked.fasta"
        applymask.main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges)
        assert_fasta(self, output_fasta, self.maskedsequence)

    def test_main_ranges_nozip(self):
        fasta_filepath = "data/test.fasta"
        mask_filepath = "data/mask_range.tsv"
        mask_format = "range"
        use_gzip = "false"
        print_mask_ranges = "true"
        output_fasta = "test.masked.fasta"
        applymask.main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges)
        assert_fasta(self, output_fasta, self.maskedsequence)


    def test_main_fasta_zip(self):
        fasta_filepath = "data/test.fasta.gz"
        mask_filepath = "data/mask_fasta.fasta"
        mask_format = "fasta"
        use_gzip = "true"
        print_mask_ranges = "true"
        output_fasta = "test.fasta.masked.gz"
        applymask.main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges)
        assert_fasta_zip(self, output_fasta, self.maskedsequence)

    def test_main_position_zip(self):
        fasta_filepath = "data/test.fasta.gz"
        mask_filepath = "data/mask_position.txt"
        mask_format = "position"
        use_gzip = "true"
        print_mask_ranges = "true"
        output_fasta = "test.fasta.masked.gz"
        applymask.main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges)
        assert_fasta_zip(self, output_fasta, self.maskedsequence)


    def test_main_ranges_zip(self):
        fasta_filepath = "data/test.fasta.gz"
        mask_filepath = "data/mask_range.tsv"
        mask_format = "range"
        use_gzip = "true"
        print_mask_ranges = "true"
        output_fasta = "test.fasta.masked.gz"
        applymask.main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges)
        assert_fasta_zip(self, output_fasta, self.maskedsequence)

    def test_main_unknown_mask(self):
        fasta_filepath = "data/test.fasta.gz"
        mask_filepath = "data/mask_range.tsv"
        mask_format = "something"
        use_gzip = "true"
        print_mask_ranges = "true"
        output_fasta = "test.fasta.masked.gz"
        applymask.main(mask_filepath, mask_format, fasta_filepath, use_gzip, print_mask_ranges)
        assert_fasta_zip(self, output_fasta, self.maskedsequence)



if __name__ == "__main__":
    unittest.main()