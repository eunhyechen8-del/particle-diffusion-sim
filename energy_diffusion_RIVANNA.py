import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import random
job_id = os.getenv("SLURM_JOB_ID", "local")
#notes from prof: shorten code by removing repetitive code, especiailly with the case handling
# get rivanna working, get in contact with the grad student? also i can still have trial limit instead of infinite, just make it a super big number lol
# #energy_diff = energy_cost - ( V * (neighbor_count(particle) + neighbor_count((row,(column - 1) % L)) - 1)) this part of the annailation is unused BUT technically redundant bc if u apply acceptnace prob it WILL ALWAYS ACCEPT bc its more favorable due to the acceptance prob

V = 1
energy_cost = -10
T = 20
L = 100 #size of shape sides #12 FOR THE MOVIE
trials = 500      # number of trials to run on my laptop
t = 2000      # Time steps per trial, obviously increase time if L is greater #200 for MOVIE
initial_prob = 0.3
N = L**2
initialAmount = int(initial_prob * N)  # +1 to make it even for the movie
print(initialAmount)

def neighbor_count(site,currentLattice):
    row_site = site[0]
    column_site = site[1]
    num_neigbors = 0
    if currentLattice[row_site][(column_site - 1) % L] == 1:
        num_neigbors += 1
    if currentLattice[row_site][(column_site + 1) % L] == 1:
        num_neigbors += 1
    if currentLattice[(row_site + 1) % L][column_site] == 1:
        num_neigbors += 1
    if currentLattice[(row_site - 1) % L][column_site] == 1:
        num_neigbors += 1
    if currentLattice[(row_site - 1) % L][(column_site - 1) % L] == 1:
        num_neigbors += 1
    if currentLattice[(row_site + 1) % L][(column_site + 1) % L] == 1:
        num_neigbors += 1
    return num_neigbors

def acceptance(energy):
    A =  min(1,math.e**(-energy/T))
    rand_num = random.uniform(0,1)
    if rand_num < A:
        return True
    return False

def total_energy_configuration(initial_lattice):
    count = 0
    particles_in_config = np.argwhere(initial_lattice == 1)
    for site in particles_in_config:
        count += neighbor_count(site, initial_lattice)
    return V * (count/2)

all_trials_results = [] #collects every single trial history, we want this to take the average <N> vs t, so we can get the smooth curve for 1/t
all_energy_history = []
average_energy_history = []
all_energy_per_particle = []
for trial in range(trials):
    #resets
    trial_energy_per_particle= [] #used to graph E/N later, shows how the on aervage particles have less and less neighbors, stays true if annahlation cost is 0
    total_energy_per_step = [] #records the total energy history of a single trial at every time step, temporary and is refreshed every trial,
    total_energy = 0
    energy_diff = 0
    lattice = np.zeros((L, L), dtype=int)
    particle_count = 0
    each_trial_history = []
    plotTrial = [] # this is where im keeping my data

    #throwing particles into LxL neighborhood config
    while particle_count < initialAmount:
        row = random.randint(0,L-1)
        column = random.randint(0,L-1)
        if lattice[row][column] == 0:
            lattice[row][column] = 1
            particle_count += 1

    #THESE TWO LINES TRACKS THE INITIAL CONDITIONS , may have been the reason why my program ended showing two ghost particles
    each_trial_history.append(int(np.sum(lattice)))  # keeps track of HOW MANY particles at each time step
    plotTrial.append(lattice.copy())  # saves a copy of the current lattice at this time step into the list, for caluclating position of every particle
    total_energy = total_energy_configuration(lattice)
    total_energy_per_step.append(total_energy)
    trial_energy_per_particle.append(total_energy_configuration(lattice) /np.sum(lattice))

    for time_step in range(t):
        current_positions = np.argwhere(lattice)

        np.random.shuffle(current_positions)

        for particle in current_positions:
            energy_diff = 0
            r = random.randint(0, 6)
            row = particle[0]
            column = particle[1]
            new_pos = ()
            if lattice[row][
                column] == 0:  # checks if particle is still there because theres a possiblity that the particle may have been annailated by a previous particle
                continue
            match r:
                case 0:  # does not move
                    continue
                case 1:  # moves to the left. so column - 1
                    new_pos = (row, (column - 1) % L)
                case 2:  # moves to the right, so column + 1
                    new_pos = (row, (column + 1) % L)
                case 3:  # moves down so row + 1
                    new_pos = ((row + 1) % L, column)
                case 4:  # moves up so row - 1
                    new_pos = ((row - 1) % L, column)
                case 5:  # moves up and left so row and column both -1
                    new_pos = ((row - 1) % L, (column - 1) % L)
                case 6:  # moves down and right so row and column both +1
                    new_pos = ((row + 1) % L, (column + 1) % L)

            new_row = new_pos[0]
            new_column = new_pos[1]
            if lattice[new_row][new_column] == 1:
                energy_diff = energy_cost + V * ((neighbor_count((new_row, new_column), lattice)) + neighbor_count(particle, lattice) - 1) # all of the neighbors of BOTH particles will lose bond with the annahihlated particles, and we need to account for the fact both will count the shared bond between the two particles
                if acceptance(energy_diff):  # technically this will ALWAYS accept, which is what we want
                    lattice[row][column] = 0
                    lattice[new_row][new_column] = 0
                else:
                    continue  # move was not accepted, so particle does not move, technically unreachable in this case bc annaihlation is always favorable
            else:
                energy_diff = V * ((neighbor_count((new_row, new_column), lattice)) - neighbor_count(particle, lattice))
                if acceptance(energy_diff):
                    lattice[new_row][new_column] = 1
                    lattice[row][column] = 0
                else:
                    continue  # move not accepted so we stay in the same spot

        current_energy = total_energy_configuration(lattice) #gets the total energy after each step
        total_energy_per_step.append(current_energy)
        current_num_particles = np.sum(lattice)
        if current_num_particles > 0:
            trial_energy_per_particle.append(current_energy/current_num_particles)
        else:
            trial_energy_per_particle.append(0)
        each_trial_history.append(int(np.sum(lattice)))  # keeps track of HOW MANY particles at each time step
        plotTrial.append(lattice.copy())  # saves a copy of the current lattice at this time step into the list, for caluclating position of every particle
        #total_energy_per_step.append(total_energy) #records the total energy at this time step after all particles has moved
        #average_energy_history.append(total_energy_per_step / np.sum(lattice))

    all_trials_results.append(each_trial_history)
    all_energy_history.append(total_energy_per_step) #saves the history of all trials so we can average them after
    all_energy_per_particle.append(trial_energy_per_particle)
    #np.save("acceptanceMovie.npy", np.array(plotTrial)) #saving the data file

final_data = np.array(all_energy_history)
time_axis = np.arange(t+1) # includes the initial condition

# Staircase
plt.figure(figsize=(10, 5))
for i in range(min(10, trials)):
    plt.step(time_axis, final_data[i], where='post', alpha=0.5, label=f'Trial {i+1}')
plt.title("Individual Trials: Total Energy of Configuration")
plt.xlabel("Time")
plt.ylabel("Total E")
plt.xlim(0, 100)
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig(f"T=20_total_energy_{job_id}.png")
plt.close()


# Average
plt.figure(figsize=(10, 5))
plt.plot(time_axis, np.abs(np.mean(all_energy_history, axis=0)), color='black', label="Average E")
plt.plot(time_axis[1:],(1/(0.01*time_axis))[1:], label="1/t")
plt.title("Average Energy Diffusion Curve")
plt.xscale('log')   # make x-axis logarithmic
plt.yscale('log')   # make y-axis logarithmic
plt.xlabel("Time")
plt.ylabel("Avg E")
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig(f"T=20_average_total_energy{job_id}.png")
plt.close()

final_data = np.array(all_energy_per_particle)
time_axis = np.arange(t+1) # includes the initial condition

# i think this shows that the particles are moving away from each other
plt.figure(figsize=(10, 5))
plt.plot(time_axis, np.mean(final_data, axis=0), color='black', label="Average E/N")
plt.plot(time_axis, np.mean(all_energy_per_particle,axis=0)[0]/(time_axis + 1), color = "red" , label="Theory (1/t)")
plt.title("Average Energy Per Particle")
plt.xlabel("Time")
plt.ylabel("Avg E Per Particle")
plt.xlim(0, 150)
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig(f"T=20_avg_energy_per_particle_linear_{job_id}.png")
plt.close()

plt.figure(figsize=(10, 5))
plt.plot(time_axis, np.mean(final_data, axis=0), color='black', label="Average E/N")
plt.plot(time_axis, np.mean(all_energy_per_particle,axis=0)[0]/(time_axis + 1), color = "red" , label="Theory (1/t)")
plt.title("Average Energy Per Particle (LOG)")
plt.xlabel("Time")
plt.ylabel("Avg E Per Particle")
plt.legend()
plt.xscale('log')   # make x-axis logarithmic
plt.yscale('log')   # make y-axis logarithmic
plt.grid(True)
#plt.show()
plt.savefig(f"T=20_avg_energy_per_particle_loglog_{job_id}.png")
plt.close()


final_data = np.array(all_trials_results)
time_axis = np.arange(t+1) # includes the initial condition

# Staircase
plt.figure(figsize=(10, 5))
for i in range(min(10, trials)):
    plt.step(time_axis, final_data[i], where='post', alpha=0.5, label=f'Trial {i+1}')
plt.title("Individual Trials: The Staircase Pattern")
plt.xlabel("Time")
plt.ylabel("N")
plt.xlim(0, 200)
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig(f"T=20_particle_count_staircase_{job_id}.png")
plt.close()

# average curve <N>
plt.figure(figsize=(10, 6))
avg_n = np.mean(final_data, axis=0)
plt.plot(time_axis, avg_n, color='black', linewidth=2)
plt.title(f"Average Particle Decay over {len(final_data)} Trials")
plt.xlabel("Time (t)")
plt.ylabel("Average N")
plt.grid(True, alpha=0.2)
# plt.show()
plt.savefig(f"T=20_average_particle_count_{job_id}.png")
plt.close()

plt.figure(figsize=(10, 5))
plt.plot(time_axis, np.mean(final_data, axis=0), color='black', label="Average N")
plt.plot(time_axis,(1/(0.04*time_axis)), label="1/t")
plt.xscale('log')   # make x-axis logarithmic
plt.yscale('log')   # make y-axis logarithmic
plt.title("Average Number of Particles Decay LOG Curve")
plt.xlabel("Time")
plt.ylabel("Avg N")
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig(f"T=20_LOG_average_particle_count_{job_id}.png")
plt.close()
