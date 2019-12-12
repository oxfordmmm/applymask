# Test applymask.py
# Run all tests: python3 test_applymask.py
# Run one test:  python3 test_applymask.py TestApplyMask.test_load_mask_format1
# Only skip tests for debugging reason. All tests should be passed

import io

import unittest
import applymask

class TestApplyMask(unittest.TestCase):
    def setUp(self):
        self.fasta = "tests/test.fasta"
        self.use_gzip = "false"
        self.print_mask_ranges = "true"

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

    def test_load_mask_range(self):
        data = ["0\t0","5\t7", "9\t9"]
        expected = list("1000011101")

        result = applymask.load_mask_range_aux(data, 10)
        self.assertEqual(result, expected)

    def test_load_mask_position_aux(self):
        data = ["0","5", "6", "7", "9"]
        expected = list("1000011101")

        result = applymask.load_mask_position_aux(data, 10)
        self.assertEqual(result, expected)
    
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

if __name__ == "__main__":
    unittest.main()