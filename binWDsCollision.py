from loadmodules import *

name1 = "../output08/snapshot08_010"
name2 = "../output08/snapshot08_010"

binname = "bin.dat.ic"
rhocut = 1.0 # g/ccm
velocity = 8e8 #8000 km/s
overlap_factor = 2.0/3.0

sp = loaders.load_species("species55.txt")
iC12 = sp['names'].index('c12')
iO16 = sp['names'].index('o16')


def compute_cumulative_mass(s):
    rsort = s.r().argsort()

    mcum = np.zeros(s.npart)
    mcum[0] = s.mass[rsort[0]]
    for i in range(1, s.npart):
        mcum[rsort[i]] = mcum[rsort[i - 1]] + s.mass[rsort[i]]
    s.data['mcum'] = mcum
    return


s1 = gadget_readsnapname(name1, hdf5=True, loadonlytype=[0])
s2 = gadget_readsnapname(name2, hdf5=True, loadonlytype=[0])

compute_cumulative_mass(s1)
compute_cumulative_mass(s2)

i1, = np.where(s1.rho > rhocut)
i2, = np.where(s2.rho > rhocut)

m1 = s1.data["mass"][i1].astype('float64').sum()
m2 = s2.data["mass"][i2].astype('float64').sum()
m = m1 + m2

r1 = s1.r()[i1].max()
r2 = s2.r()[i2].max()
print("R1:", r1)
print("R2:", r2)


q = m2 / m1

b = (r1 + r2) * (1 - overlap_factor)

if r1 > r2:
    a0 = r1 * (0.6 * q ** (2. / 3.) + log(1 + q ** (1. / 3.))) / (0.49 * q ** (2. / 3.))
else:
    a0 = r2 * (0.6 * q ** (2. / 3.) + log(1 + q ** (1. / 3.))) / (0.49 * q ** (2. / 3.))

x1 = m2 / m * a0
x2 = m1 / m * a0

y1 = m2 / m * b
y2 = m1 / m * b

voffset1 = zeros((3, s1.npart))
voffset2 = zeros((3, s2.npart))

npart1 = size(i1)
npart2 = size(i2)
npart = npart1 + npart2

data = {}
data['pos'] = np.zeros((npart, 3))
data['vel'] = np.zeros((npart, 3))
data['bfld'] = np.zeros((npart, 3))
data['mass'] = np.zeros(npart)
data['u'] = np.zeros(npart)
data['xnuc'] = np.zeros((npart, sp['count']))
data['pass'] = np.zeros((npart, 2))
data['count'] = npart

c1 = (s1.pos[i1, :].astype('f8') * s1.data['mass'][i1][:, None]).sum(axis=0) / m1
c1 = -c1
c1[0] -= x1
c1[1] -= y1

data['pos'][:npart1, :] = s1.pos[i1, :] + c1[None, :]
data['vel'][:npart1, :] = s1.vel[i1, :]
data['mass'][:npart1] = s1.mass[i1]
data['u'][:npart1] = s1.data['u'][i1]
data['xnuc'][:npart1, :] = s1.data['xnuc'][i1, :]
data['pass'][:npart1, 0] = 1.0

c2 = (s2.pos[i2, :].astype('f8') * s2.data['mass'][i2][:, None]).sum(axis=0) / m2
c2 = -c2
c2[0] += x2
c2[1] += y2

data['pos'][npart1:, :] = s2.pos[i2, :] + c2[None, :]
data['vel'][npart1:, :] = s2.vel[i2, :]
data['mass'][npart1:] = s2.mass[i2]
data['u'][npart1:] = s2.data['u'][i2]
data['xnuc'][npart1:, :] = s2.data['xnuc'][i2, :]
data['pass'][npart1:, 1] = 1.0

data['vel'][:npart1, 0] = data['vel'][:npart1, 0] + velocity #some velocity in the x direction
data['vel'][npart1:, 0] = data['vel'][npart1:, 0] - velocity

mm = np.array([0., 0., 1e3 * 1e9 ** 3 / 2.])  # 1e3 G at 1e9 cm

c1 = (data['pos'] * data['mass'][:, None] * data['pass'][:, 0][:, None]).sum(axis=0) / m1
r1 = np.maximum(np.sqrt(((data['pos'] - c1[None, :]) ** 2).sum(axis=1)), 3e7).astype('float64')
i, = np.where(r1 < 1e10)
rad1 = data['pos'][i, :] - c1[None, :]
data['bfld'][i, :] = 3. * rad1 * (mm[None, :] * rad1).sum(axis=1)[:, None] / (r1[i] ** 5)[:, None] \
                     - mm[None, :] / (r1[i] ** 3)[:,None]

c2 = (data['pos'] * data['mass'][:, None] * data['pass'][:, 1][:, None]).sum(axis=0) / m2
r2 = np.maximum(np.sqrt(((data['pos'] - c2[None, :]) ** 2).sum(axis=1)), 3e7).astype('float64')
i, = np.where(r2 < 1e10)
rad2 = data['pos'][i, :] - c2[None, :]
data['bfld'][i, :] += 3. * rad2 * (mm[None, :] * rad2).sum(axis=1)[:, None] / (r2[i] ** 5)[:, None] \
                    - mm[None, :] / (r2[i] ** 3)[:,None]

#TODO: Check many different uses of r1 and r2!!

print(c1)
print(c2)

xnuc = np.zeros(sp['count'])
xnuc[iC12] = 0.5
xnuc[iO16] = 0.5

data['boxsize'] = 1e10
data['pos'] += 0.5 * data['boxsize']
gadget_add_grids(data, [1e10, 1e11, 1e12], 32, xnuc=xnuc)

gadget_write_ics(binname, data, double=True, format="hdf5")

s = gadget_readsnapname(binname)
s.plot_pos()