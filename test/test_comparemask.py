# Test comparemask.py
# Run all tests: python3 test_comparemask.py
# Run code coverage: coverge run test_comparemask.py 
# View code coverage report: coverage report -m
# Generate code coverage html report: coverage html

import os
import unittest
import comparemask

def assert_mask(self, mask_name):
    mask = comparemask.load_position_mask(mask_name)
    return mask

class TestApplyMask(unittest.TestCase):
    def setUp(self):
        self.expected_mask_positions = [0,59,60,119,120,179,180,239,240,299]
        self.expected_mask_strings = '0\n59\n60\n119\n120\n179\n180\n239\n240\n299'
        self.mask1 = [1,3,5,7,9]
        self.mask2 = [2,3,5,7,11]
        
    def test_load_position_mask(self):
        expected = self.expected_mask_positions
        result = assert_mask(self,"tests/mask_position.txt")
        self.assertEqual(expected, result)
    
    def test_compare_position_mask(self):
        expected_diff_1_2 = [1,9]
        expected_diff_2_1 = [2,11]
        expected_diff_all = [1,2,9,11]
        result1, result2, result_all = comparemask.compare_position_mask(self.mask1, self.mask2)
        self.assertEqual(expected_diff_1_2, result1)
        self.assertEqual(expected_diff_2_1, result2)
        self.assertEqual(expected_diff_all, result_all)

    def test_fasta_to_posistions(self):
        expected = self.expected_mask_positions
        result = comparemask.fasta_to_posistions("tests/mask_fasta.fasta")
        self.assertEqual(expected, result)

    def test_write_position_to_file(self):
        expected = self.expected_mask_strings
        result = comparemask.write_position_to_file(self.expected_mask_positions, 'tests/output_positions.txt')
        self.assertEqual(expected, result)


if __name__ == "__main__":
    unittest.main()