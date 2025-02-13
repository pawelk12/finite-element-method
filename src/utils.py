import sys

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
