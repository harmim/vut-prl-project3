/**
 * VUT FIT PRL 2020 Project - Line-of-Sight.
 *
 * @author Dominik Harmim <harmim6@gmail.com>
 */


#include <cstdlib>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <limits>
#include <climits>
#include <mpi.h>

#ifdef DEBUG
#include <chrono>
#endif


#define COMM MPI_COMM_WORLD /// Default MPI communicator.
#define MASTER 0 /// A rank of the master process.
#define TAG 0 /// An MPI tag used for the transmission of messages.
#define EMPTY INT_MIN /// An empty value of an MPI message.


using namespace std;


/**
 * Prints an error message due to an MPI error and terminates the program with
 * an erroneous exit code.
 */
auto MPI_error() -> void
{
	cerr << "Error: an MPI library call has failed." << endl;
	MPI_Abort(COMM, EXIT_FAILURE);
}


/**
 * The master process parses and validates input altitudes, checks the current
 * number of processes, and sends the altitudes to the others processes.
 *
 * @param rank A rank of a process.
 * @param procs_count The number of all processes.
 * @param argc The number of input arguments.
 * @param argv An array of input arguments.
 * @return Returns the number of altitudes (expect the first one - a place of
 *         observation).
 */
auto send_alts(
	const int rank,
	const int procs_count,
	const int argc,
	const char *const argv[]
) -> size_t
{
	if (rank != MASTER)
	{
		return 0;
	}

	if (argc != 2)
	{
		cerr << "Error: expecting one argument: " << argv[0] << " altitudes"
			<< endl;
		MPI_Abort(COMM, EXIT_FAILURE);
	}

	if (!regex_match(argv[1], regex("^-?[0-9]+(,-?[0-9]+)*$")))
	{
		cerr << "Error: invalid format of altitudes,"
			" expecting '^-?[0-9]+(,-?[0-9]+)*$'." << endl;
		MPI_Abort(COMM, EXIT_FAILURE);
	}

	vector<int> alts;
	stringstream alts_stream(argv[1]);
	for (string alt; getline(alts_stream, alt, ','); )
	{
		alts.push_back(stoi(alt));
	}
	const size_t alts_count = alts.size() - 1;

	int expected_procs_count = alts_count;
	if (alts_count > 1)
	{
		expected_procs_count =
			static_cast<int>(pow(2, ceil(log(alts_count) / log(2))) / 2);
	}
	if (alts_count > 0 && procs_count != expected_procs_count)
	{
		cerr << "Error: the expected number of processes is "
			<< expected_procs_count << " but the current number of processes"
			" is " << procs_count << "." << endl;
		MPI_Abort(COMM, EXIT_FAILURE);
	}

	if (!alts_count)
	{
		return alts_count;
	}

	auto alts_iterator = ++alts.begin();
	for (int r = MASTER; r < procs_count; r++)
	{
		for (size_t i = 0; i < 3; i++)
		{
			int alt = EMPTY;
			if (!i)
			{
				alt = alts[0];
			}
			else if (alts_iterator != alts.end())
			{
				alt = *alts_iterator++;
			}
			if (MPI_Send(&alt, 1, MPI_INT, r, TAG, COMM))
			{
				MPI_error();
			}
		}
	}

	return alts_count;
}


/**
 * Receives altitudes from the master process.
 *
 * @param alts Received altitudes.
 */
auto receive_alts(int alts[]) -> void
{
	for (size_t i = 0; i < 3; i++)
	{
		if (MPI_Recv(&alts[i], 1, MPI_INT, MASTER, TAG, COMM, nullptr))
		{
			MPI_error();
		}
	}
}


/**
 * Calculates vertical angles from altitudes.
 *
 * @param angles Calculated vertical angles.
 * @param alts Altitudes for the calculation of vertical angles.
 * @param rank A rank of a process.
 */
auto calculate_angles(double angles[], const int alts[], const int rank) -> void
{
	for (size_t i = 0; i < 2; i++)
	{
		angles[i] = numeric_limits<double>::lowest();
		if (alts[i + 1] != EMPTY)
		{
			const auto length = static_cast<double>(rank * 2 + 1 + i);
			angles[i] = atan((alts[i + 1] - alts[0]) / length);
		}
	}
}


/**
 * Calculates maximal previous vertical angles (prescan).
 *
 * @param max_angles Calculated maximal previous vertical angles.
 * @param angles Vertical angles of points.
 * @param rank A rank of a process.
 * @param procs_count The number of all processes.
 */
auto max_prescan(
	double max_angles[],
	const double angles[],
	const int rank,
	const int procs_count
) -> void
{
	max_angles[0] = angles[0];
	max_angles[1] = angles[1];

	if (procs_count == 1)
	{
		max_angles[1] = max_angles[0];
		max_angles[0] = numeric_limits<double>::lowest();
		return;
	}

	int step, prev_step, r;
	double neigh;

	// up-sweep
	max_angles[1] = max(max_angles[0], max_angles[1]);
	for (int d = 0; d <= log2(procs_count) - 1; d++)
	{
		step = pow(2, d + 1);
		prev_step = pow(2, d);

		if (!((rank + 1) % step))
		{
			r = rank - prev_step;
			if (MPI_Recv(&neigh, 1, MPI_DOUBLE, r, TAG, COMM, nullptr))
			{
				MPI_error();
			}
			max_angles[1] = max(max_angles[1], neigh);
		}
		else if (!((rank + 1) % prev_step))
		{
			r = rank + prev_step;
			if (MPI_Send(&max_angles[1], 1, MPI_DOUBLE, r, TAG, COMM))
			{
				MPI_error();
			}
		}
	}

	// down-sweep
	if (rank == procs_count - 1)
	{
		max_angles[1] = numeric_limits<double>::lowest();
	}
	for (int d = static_cast<int>(log2(procs_count)) - 1; d >= 0; d--)
	{
		step = pow(2, d + 1);
		prev_step = pow(2, d);

		if (!((rank + 1) % step))
		{
			r = rank - prev_step;
			if (MPI_Send(&max_angles[1], 1, MPI_DOUBLE, r, TAG, COMM))
			{
				MPI_error();
			}
			if (MPI_Recv(&neigh, 1, MPI_DOUBLE, r, TAG, COMM, nullptr))
			{
				MPI_error();
			}
			max_angles[1] = max(max_angles[1], neigh);
		}
		else if (!((rank + 1) % prev_step))
		{
			r = rank + prev_step;
			if (MPI_Recv(&neigh, 1, MPI_DOUBLE, r, TAG, COMM, nullptr))
			{
				MPI_error();
			}
			if (MPI_Send(&max_angles[1], 1, MPI_DOUBLE, r, TAG, COMM))
			{
				MPI_error();
			}
			max_angles[1] = neigh;
		}
	}
	const double t = max_angles[0];
	max_angles[0] = max_angles[1];
	max_angles[1] = max(max_angles[1], t);
}


/**
 * Determines whether points are visible and sends this visibility to the
 * master process.
 *
 * @param angles Vertical angles of points.
 * @param max_angles Maximal previous vertical angles.
 * @param alts Altitudes of points.
 */
auto send_visibility(
	const double angles[], const double max_angles[], const int alts[]
) -> void
{
	if (alts[1] == EMPTY)
	{
		return;
	}

	for (size_t i = 0; i < 2; i++)
	{
		int visible = EMPTY;
		if (alts[i + 1] != EMPTY)
		{
			visible = angles[i] > max_angles[i];
		}
		if (MPI_Send(&visible, 1, MPI_INT, MASTER, TAG, COMM))
		{
			MPI_error();
		}
	}
}


/**
 * Receives the visibility of points from all processes.
 *
 * @param visibility The received visibility of all points.
 * @param rank A rank of a process.
 * @param alts_count The number of altitudes.
 */
auto receive_visibility(
	bool visibility[], const int rank, const size_t alts_count
) -> void
{
	if (rank != MASTER)
	{
		return;
	}

	for (size_t c = 0; c < ceil(alts_count / 2.0); c++)
	{
		for (size_t i = 0; i < 2; i++)
		{
			int visible;
			if (MPI_Recv(&visible, 1, MPI_INT, c, TAG, COMM, nullptr))
			{
				MPI_error();
			}

			if (visible == EMPTY)
			{
				break;
			}

			visibility[c * 2 + i] = visible;
		}
	}
}


/**
 * Prints given visibilities to the standard output.
 *
 * @param visibility Visibilities to be printed.
 * @param rank A rank of a process.
 * @param alts_count The number of altitudes.
 */
auto print_visibility(
	const bool visibility[], const int rank, const size_t alts_count
) -> void
{
	if (rank != MASTER)
	{
		return;
	}

	cout << "_";
	for (size_t c = 0; c < alts_count; c++)
	{
		cout << "," << (visibility[c] ? "v" : "u");
	}
	cout << endl;
}


/**
 * The entry point of the program.
 *
 * @param argc The number of input arguments.
 * @param argv An array of input arguments.
 * @return Returns EXIT_SUCCESS if the program succeeded, EXIT_FAILURE
 *         otherwise.
 */
auto main(int argc, char *argv[]) -> int
{
	if (MPI_Init(&argc, &argv))
	{
		MPI_error();
	}

	int rank, procs_count;
	if (MPI_Comm_rank(COMM, &rank))
	{
		MPI_error();
	}
	if (MPI_Comm_size(COMM, &procs_count))
	{
		MPI_error();
	}

	const auto alts_count = send_alts(rank, procs_count, argc, argv);
	if (!alts_count && rank == MASTER)
	{
		cout << "_" << endl;
		MPI_Abort(COMM, EXIT_SUCCESS);
	}

	int alts[3];
	receive_alts(alts);

#ifdef DEBUG
	chrono::time_point<chrono::high_resolution_clock> start, end;
	MPI_Barrier(COMM);
	if (rank == MASTER) start = chrono::high_resolution_clock::now();
#endif

	double angles[2];
	calculate_angles(angles, alts, rank);

	double max_angles[2];
	max_prescan(max_angles, angles, rank, procs_count);

	send_visibility(angles, max_angles, alts);

#ifdef DEBUG
	MPI_Barrier(COMM);
	if (rank == MASTER) end = chrono::high_resolution_clock::now();
#endif

	bool *const visibility = new bool[alts_count];
	receive_visibility(visibility, rank, alts_count);
	print_visibility(visibility, rank, alts_count);

#ifdef DEBUG
	if (rank == MASTER)
	{
		chrono::duration<double> diff = end - start;
		cout << "Time: " << diff.count() << endl;
	}
#endif

	delete[] visibility;
	if (MPI_Finalize())
	{
		MPI_error();
	}

	return EXIT_SUCCESS;
}
