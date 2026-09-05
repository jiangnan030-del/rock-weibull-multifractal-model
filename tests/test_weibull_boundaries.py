import unittest

from rock_weibull_model import default_demo_peak, solve_weibull_parameters


class WeibullBoundaryTests(unittest.TestCase):
    def test_lambda_one_is_valid_when_r_is_valid(self):
        result = solve_weibull_parameters(default_demo_peak(), 1.0)
        self.assertGreater(result.m, 0.0)
        self.assertGreater(result.F0, 0.0)

    def test_lambda_zero_and_above_one_are_rejected(self):
        for value in (0.0, -0.1, 1.0001):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    solve_weibull_parameters(default_demo_peak(), value)


if __name__ == "__main__":
    unittest.main()
