import sys
import math

def get_args():
    if len(sys.argv) != 3:
        print('python main.py <file_name.txt> <number_of_integration_points(2 or 3 or 4)>')
        sys.exit(1)

    filename = sys.argv[1]
    try:
        n_integration_points = int(sys.argv[2])
    except:
        print('Integration points were not provided correctly')
        sys.exit(1)

    if n_integration_points < 2 or n_integration_points > 4:
        print('Invalid number of integration points, please enter 2 or 3 or 4')
        sys.exit(1)

    return filename, n_integration_points


def init_integration_points_and_weights(n_integration_points):
    match n_integration_points:
        case 2:
            # two-point Gaussian quadrature
            integration_points = [(-1 / math.sqrt(3), -1 / math.sqrt(3)),
                                (1 / math.sqrt(3), -1 / math.sqrt(3)),
                                (1 / math.sqrt(3), 1 / math.sqrt(3)),
                                (-1 / math.sqrt(3), 1 / math.sqrt(3))]
            weights_for_integration_points = [(1, 1), (1, 1), (1, 1), (1, 1)]
            return integration_points, weights_for_integration_points

        case 3:
            # three-point Gaussian quadrature
            integration_points = [(-math.sqrt(3) / math.sqrt(5), -math.sqrt(3) / math.sqrt(5)),
                                (-math.sqrt(3) / math.sqrt(5), 0),
                                (-math.sqrt(3) / math.sqrt(5), math.sqrt(3) / math.sqrt(5)),
                                (0, -math.sqrt(3) / math.sqrt(5)),
                                (0, 0),
                                (0, math.sqrt(3) / math.sqrt(5)),
                                (math.sqrt(3) / math.sqrt(5), -math.sqrt(3) / math.sqrt(5)),
                                (math.sqrt(3) / math.sqrt(5), 0),
                                (math.sqrt(3) / math.sqrt(5), math.sqrt(3) / math.sqrt(5))]

            weights_for_integration_points = [(5 / 9, 5 / 9), (5 / 9, 8 / 9), (5 / 9, 5 / 9),
                                            (8 / 9, 5 / 9), (8 / 9, 8 / 9), (8 / 9, 5 / 9),
                                            (5 / 9, 5 / 9), (5 / 9, 8 / 9), (5 / 9, 5 / 9)]
            return integration_points, weights_for_integration_points

        case 4:
            # four-point Gaussian quadrature
            integration_points = [(-0.861136, -0.861136),
                                  (-0.861136, -0.339981),
                                  (-0.861136, 0.339981),
                                  (-0.861136, 0.861136),
                                  (-0.339981, -0.861136),
                                  (-0.339981, -0.339981),
                                  (-0.339981, 0.339981),
                                  (-0.339981, 0.861136),
                                  (0.339981, -0.861136),
                                  (0.339981, -0.339981),
                                  (0.339981, 0.339981),
                                  (0.339981, 0.861136),
                                  (0.861136, -0.861136),
                                  (0.861136, -0.339981),
                                  (0.861136, 0.339981),
                                  (0.861136, 0.861136)]

            weights_for_integration_points = [(0.347855, 0.347855), (0.347855, 0.652145), (0.347855, 0.652145),
                                              (0.347855, 0.347855),
                                              (0.652145, 0.347855), (0.652145, 0.652145), (0.652145, 0.652145),
                                              (0.652145, 0.347855),
                                              (0.652145, 0.347855), (0.652145, 0.652145), (0.652145, 0.652145),
                                              (0.652145, 0.347855),
                                              (0.347855, 0.347855), (0.347855, 0.652145), (0.347855, 0.652145),
                                              (0.347855, 0.347855)]
            return integration_points, weights_for_integration_points