'''
Created on 2014-05-01

@author: cbrunet
'''
import unittest

from fibermodes.mode import Mode, Family


class TestMode(unittest.TestCase):

    def testEquality(self):
        self.assertEqual(Mode('HE', 1, 1), Mode('HE', 1, 1))
        self.assertNotEqual(Mode('HE', 1, 1), Mode('HE', 2, 1))
        self.assertNotEqual(Mode('HE', 1, 1), Mode('LP', 0, 1))

    def testSortModes(self):
        modes = [
            Mode(Family.LP, 0, 1),
            Mode(Family.HE, 1, 1),
            Mode(Family.LP, 1, 1),
            Mode(Family.TE, 0, 1),
            Mode(Family.HE, 2, 1),
            Mode(Family.TM, 0, 1),
            Mode(Family.LP, 2, 1),
            Mode(Family.EH, 1, 1),
            Mode(Family.HE, 3, 1),
            Mode(Family.EH, 2, 1),
            Mode(Family.HE, 4, 1),
            Mode(Family.LP, 0, 2),
            Mode(Family.HE, 1, 2),
            Mode(Family.TE, 0, 2),
            Mode(Family.HE, 2, 2),
            Mode(Family.TM, 0, 2),
            Mode(Family.EH, 1, 2),
            Mode(Family.HE, 3, 2),
            Mode(Family.EH, 2, 2),
            Mode(Family.HE, 4, 2)
            ]
        for i in range(len(modes)):
            for j in range(len(modes)):
                if i < j:
                    self.assertLess(modes[i], modes[j])
                    self.assertLessEqual(modes[i], modes[j])
                    self.assertGreater(modes[j], modes[i])
                    self.assertGreaterEqual(modes[j], modes[i])
                    self.assertNotEqual(modes[i], modes[j])
                elif i > j:
                    self.assertLess(modes[j], modes[i])
                    self.assertLessEqual(modes[j], modes[i])
                    self.assertGreater(modes[i], modes[j])
                    self.assertGreaterEqual(modes[i], modes[j])
                    self.assertNotEqual(modes[j], modes[i])
                else:
                    self.assertEqual(modes[i], modes[j])
                    self.assertLessEqual(modes[i], modes[j])
                    self.assertGreaterEqual(modes[i], modes[j])


if __name__ == "__main__":
    unittest.main()
