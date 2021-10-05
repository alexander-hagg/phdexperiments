import lettuce as lt
import torch
import numpy as np
import matplotlib.pyplot as plt
import imageio
import scipy.io
from datetime import datetime
import os

# Not interactive to prevent trying to spawn a window
plt.ioff()


# get GPU ID
GPUID = os.getenv('GPUID')
print(GPUID)
deviceName = 'cuda:'+str(GPUID)
print(deviceName)

startTime = datetime.now()
# device=torch.device("cpu")
device=torch.device(deviceName)

lattice = lt.Lattice(lt.D2Q9,device=device,dtype=torch.float32)

building = imageio.imread('building.png')
building = building.astype(bool)

# Here you can change the resolution and the Reynolds number.
flow=lt.Obstacle2D(300,200,10000,0.075,lattice,building.shape[0])
numSteps = 2000  # 2000
stepSize = 50    # 50

x = flow.grid
mask_np = np.zeros([flow.resolution_x,flow.resolution_y],dtype=bool)
relative_position_x = int(mask_np.shape[0]/4-building.shape[0]/2)
relative_position_y = int(mask_np.shape[1]/2-building.shape[1]/2)

mask_np[relative_position_x:relative_position_x+building.shape[0],relative_position_y:relative_position_y+building.shape[1]]=building
mask_np=mask_np.astype(bool)
mask=mask_np.tolist()
mask=torch.tensor(mask,dtype=torch.uint8)

flow.mask=mask
collision=lt.BGKSmagorinskiCollision(lattice, tau=flow.units.relaxation_parameter_lu)
streaming=lt.StandardStreaming(lattice)

simulation=lt.Simulation(flow,lattice,collision,streaming)
MaxURep = lt.flows.obstacle.MaxUReporter(lattice,flow,numSteps/3,50)
simulation.reporters.append(MaxURep)
USampleRep = lt.flows.obstacle.uSampleReporter(lattice,flow,400)
# simulation.initialize(max_num_steps=10)

getPictures = True
enstrophies = []
print("Viscosity in lattice units:", flow.units.viscosity_lu)
print("Active nodes:", simulation.no_collision_mask[simulation.no_collision_mask is False].size())

for i in range(numSteps):
    print('Step ', i*stepSize, 'MLUPS: ', simulation.step(stepSize))
    # if i > numSteps/3:
    enstr = [(flow.calcEnstrophy(simulation.f, lattice))]
    enstrophies += enstr
    fAllEnstr = open("AllEnstrophies","a")
    fAllEnstr.write("%10.10f\n" %(enstr[0]))
    fAllEnstr.close()
    u0 = (lattice.convert_to_numpy(lattice.u(simulation.f)[0])).transpose([1, 0])
    fu0 = open("u0","a")
    fu0.write("%10.10f,%10.10f,%10.10f,%10.10f\n" %(np.max(u0), np.mean(u0), np.min(u0), np.std(u0)))
    fu0.close()
    if getPictures:
        if ((i * stepSize) % 500 == 0):
            S = lattice.shear_tensor(simulation.f)
            S /= 2.0 * lattice.rho(simulation.f) * 1.0 * lattice.cs * lattice.cs
            S = lattice.einsum('ab,ab->', [S, S])

            u0 = lattice.u(simulation.f)[0].cpu().numpy()
            u1 = lattice.u(simulation.f)[1].cpu().numpy()
            grad_u0 = np.gradient(u0)
            grad_u1 = np.gradient(u1)
            dx = flow.units.convert_length_to_pu(1.0)
            vorticity = ((grad_u0[1] - grad_u1[0]) * (grad_u0[1] - grad_u1[0]))

            plt.imshow(np.transpose(np.sqrt(u0*u0+u1*u1)))
            plt.savefig('test{}.png'.format(i * stepSize))

    if (i * stepSize) % 50 == 0:
        u0 = lattice.u(simulation.f)[0].cpu().numpy()
        u1 = lattice.u(simulation.f)[1].cpu().numpy()
        u = lattice.u(simulation.f).cpu().numpy()
        field = np.transpose(np.sqrt(u0*u0+u1*u1))
        np.save('fig{}.npy'.format(i*stepSize), u)
        scipy.io.savemat('velocity{}.mat'.format(stepSize*i),{'u':field})

maxEnstrophies = np.mean(np.asarray(enstrophies))
maxU = np.mean(np.asarray(simulation.reporters[0].out))
print('Simulation successful: ')
print('Mean Enstrophy: ', maxEnstrophies)
print('Mean U_max: ', maxU)

# plt.plot(enstrophies)
# plt.savefig('enstrophies.pdf')
# np.save('enstropies.npy',np.asarray(enstrophies))

# plt.plot(np.asarray(simulation.reporters[0].out))
# plt.savefig('maxU.pdf')
# np.save('maxU.npy',np.asarray(simulation.reporters[0].out))

fuMax = open("uMax","a")
fenstr = open("enstrophy","a")
fuMax.write("%10.10f\n" %(maxU))
fuMax.close()
fenstr.write("%10.10f\n" %(maxEnstrophies))
fenstr.close()

# fu0.flush()
# fuMax.flush()
# fenstr.flush()
# fAllEnstr.flush()




print(datetime.now() - startTime)

sigFinished = open("done","w")
sigFinished.close()
