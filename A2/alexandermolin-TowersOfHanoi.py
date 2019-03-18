import argparse
import sys

parser=argparse.ArgumentParser(description="Tower of Hanoi solution")
parser.add_argument('-n','--number_of_disk', type = int, required=True, help = "the number of the disk")
#group = parser.add_mutually_exclusive_group()
#group.add_argument('-all', '--allnumbers', action='store_true',help="returns a comma separated list of all fibo-numbers")
args = parser.parse_args()
n = args.number_of_disk

i = 0 

def tower_of_hanoi (n, source, auxiliary, destination):
    if n == 1:
        print(f"Move Disk from {source} to {destination}")
        global i
        i +=1
    else:
        tower_of_hanoi (n-1, source, destination, auxiliary)
        print(f"Move Disk from {source} to {destination}")
        i += 1
        tower_of_hanoi (n-1, auxiliary, source, destination)
   

tower_of_hanoi(n, 'A', 'B', 'C')
print(i, file=sys.stderr)

#if args.allnumbers:
 # print(printAll(n))

#else:
 #  print(fibonacci_function(n))