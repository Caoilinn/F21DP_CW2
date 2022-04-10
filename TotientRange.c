// TotientRance.c - Sequential Euler Totient Function (C Version)
// compile: gcc -Wall -O -o TotientRange TotientRange.c
// run:     ./TotientRange lower_num uppper_num

// Greg Michaelson 14/10/2003
// Patrick Maier   29/01/2010 [enforced ANSI C compliance]

// This program calculates the sum of the totients between a lower and an
// upper limit using C longs. It is based on earlier work by:
// Phil Trinder, Nathan Charles, Hans-Wolfgang Loidl and Colin Runciman

#include <stdio.h>
#include <mpi.h>
// hcf x 0 = x
// hcf x y = hcf y (rem x y)

long hcf(long x, long y)
{
  long t;

  while (y != 0)
  {
    t = x % y;
    x = y;
    y = t;
  }
  return x;
}

// relprime x y = hcf x y == 1

int relprime(long x, long y)
{
  return hcf(x, y) == 1;
}

// euler n = length (filter (relprime n) [1 .. n-1])

long euler(long n)
{
  long length, i;

  length = 0;
  for (i = 1; i < n; i++)
    if (relprime(n, i))
      length++;
  return length;
}

// sumTotient lower upper = sum (map euler [lower, lower+1 .. upper])

long sumTotient(long lower, long upper, int mpi)
{
  MPI_Status Stat;
  int rank, world_size;
  long sum, i, part;

  // Rank(id) of current process and the overall number of all processes
  mpi = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  mpi = MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (world_size > 1)
  {
    // If the rank(process id) is 0 then it is the main process
    if (rank == 0)
    {

      // This is how many different values there are in the range
      //  i.e if the range is 5...10 then the size is 6 (10-5) + 1 = 6
      //  or 10..17 (17-10) + 1 = 8
      long size = upper - lower + 1;
      // World size - 1 becuase the 0th process isn't being used here
      long chunk_size = size / (world_size - 1);
      int remainder = size % (world_size - 1);
      int count = 1;
      // Start at 1 because you're working partial figures for each process not including the main process
      for (i = lower; i <= upper - 1; i += chunk_size)
      {
        long start, end;
        if (i == lower)
        {
          start = lower;
        }
        else
        {
          start = i;
        }

        end = start + chunk_size - 1;

        if (i == upper - chunk_size)
        {
          end = upper;
        }

        mpi = MPI_Send(&start, 1, MPI_LONG, count, 1, MPI_COMM_WORLD);
        mpi = MPI_Send(&end, 1, MPI_LONG, count, 1, MPI_COMM_WORLD);
        count = count + 1;
      }

      // This is where the sum of the totient range gets calculated
      sum = 0;

      for (i = 1; i < world_size; i++)
      {
        mpi = MPI_Recv(&part, 1, MPI_LONG, i, 2, MPI_COMM_WORLD, &Stat);
        sum += part;
      }
    }
    else
    {
      mpi = MPI_Recv(&lower, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD, &Stat);
      mpi = MPI_Recv(&upper, 1, MPI_LONG, 0, 1, MPI_COMM_WORLD, &Stat);

      part = 0;

      for (i = lower; i <= upper; i++)
      {
        part = part + euler(i);
      }
      mpi = MPI_Send(&part, 1, MPI_LONG, 0, 2, MPI_COMM_WORLD);
    }
  }
  else
  {
    for (i = lower; i <= upper; i++)
      sum = sum + euler(i);
  }
  mpi = MPI_Finalize();
  return sum;
}

int main(int argc, char **argv)
{
  long lower, upper, result = 0;

  if (argc != 3)
  {
    printf("You need to pass two arguments i.e. a lower and upper boundary\n");
    return 1;
  }

  int mpi;
  sscanf(argv[1], "%ld", &lower);
  sscanf(argv[2], "%ld", &upper);

  mpi = MPI_Init(NULL, NULL);
  result = sumTotient(lower, upper, mpi);
  // demoFunc(mpi);
  printf("Sum of Totients  between [%ld..%ld] is %ld\n",
         lower, upper, result);
  return 0;
}