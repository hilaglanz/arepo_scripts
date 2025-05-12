import loadmodules
from loadmodules import *
import argparse


def create_co_wd(mass=1.10, he4_mass=0.03, species_file="species55.txt", eos_file="/eostable/helm_table.dat", ic_file="wd.dat.ic"):
    global sp, eos, k, npart, boxsize, data, i, rho, pres, index, temp, r
    mtot = mass * msol
    mhe4 = he4_mass * msol
    sp = loaders.load_species(species_file)
    iHe4 = sp['names'].index('he4')
    iC12 = sp['names'].index('c12')
    iO16 = sp['names'].index('o16')
    eos = loadhelm_eos(ppath + eos_file, species_file, True)
    rhoc = wdCOgetRhoCFromMassExact(mtot, eos)
    wd = create_wd(eos, rhoc, xC12=0.5, xO16=0.5, tol=1e-10)
    mcum = np.cumsum(wd['dm'])
    wd_mass = wd['dm'].sum() / msol
    co_mass = wd_mass
    if mhe4 > 0:
        j, = np.where((mcum > (mcum.max() - mhe4)))
        rad_heshell = wd['r'][j[0]].min()
        k, = np.where((mcum <= (mcum.max() - mhe4)))
        print("Radius of He shell: %gkm" % (rad_heshell * 1e-5))
        print("Mass of He shell: %gMsun" % (wd['dm'][j].sum() / msol))
        co_mass = wd['dm'][k].sum() / msol
    print("Mass of WD: %gMsun" % (wd_mass))
    print("Mass of CO core: %gMsun" % (co_mass))
    sp = loaders.load_species(species_file)
    wd['v'] = np.zeros(wd['ncells'])
    wd['xnuc'] = np.zeros(sp['count'])
    wd['xnuc'][iC12] = 0.5
    wd['xnuc'][iO16] = 0.5
    wd['count'] = wd['ncells']
    npart = 0  # not set
    pmass = 1e-7 * msol
    boxsize = 1e10
    data = create_particles_healpix(wd, eos, npart, nspecies=sp['count'], minenergy=1e14, randomizeshells=True,
                                    randomizeradii=True, pmass=pmass, makebox=True, boxsize=boxsize, boxres=32,
                                    boxfactor=10.)
    frho = interpolate.interp1d(wd['r'], wd['rho'], kind='cubic')
    fpres = interpolate.interp1d(wd['r'], wd['p'], kind='cubic')
    rad = np.sqrt(((data['pos'] - 0.5 * boxsize) ** 2.).sum(axis=1))
    i, = np.where(rad < wd['r'].max())
    rho = frho(rad[i])
    pres = fpres(rad[i])
    data['xnuc'] = np.zeros((data['count'], sp['count']))
    data['xnuc'][:, iC12] = 0.5
    data['xnuc'][:, iO16] = 0.5
    if mhe4 > 0:
        j, = np.where((rad < wd['r'].max()) & (rad > rad_heshell))
        print("He shell is made up of %d cells." % size(j))
        data['xnuc'][j, iC12] = 0.
        data['xnuc'][j, iO16] = 0.
        data['xnuc'][j, iHe4] = 1.
    for index in range(size(i)):
        idx = i[index]
        temp, data['u'][idx] = eos.pgiven(rho[index], data['xnuc'][index, :], pres[index])
    data['mass'][:] = 0.
    data['mass'][i] = rho
    data['bfld'] = np.zeros((data['count'], 3))
    data['pass'] = np.zeros((data['count'], 2))
    if mhe4 > 0:
        data['pass'][j, 1] = 1.
    data['pass'][i, 0] = 1.
    r = rad[i].max()
    print("WD Radius: %gkm" % (r * 1e-5))
    gadget_write_ics(ic_file, data, double=False, format="hdf5")





def InitParser():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--mesa_file', type=str, help='path to mesa model file', default="")
    parser.add_argument('--ic_file_name', type=str, help='path to save the ic file', default="bin.ic.dat")
    parser.add_argument('--wd_mass', type=float, help='if no mesa file create a model with this mass in msol', default=1.1)
    parser.add_argument('--he4_mass', type=float, help='if no mesa file create a model with this he4 mass in msol', default=0.03)
    parser.add_argument('--species_file', type=str, help='path to species file', default="species55.txt")
    parser.add_argument('--eos_file', type=str, help='path to eos file', default="helm_table.dat")
    return parser

if __name__ == "__main__":
    for arg in sys.argv:
        print(arg)
    print(len(sys.argv))
    parser = InitParser()
    args = parser.parse_args()
    create_co_wd(mass=args.wd_mass, he4_mass=args.he4_mass, species_file=args.species_file,
                 eos_file=args.eos_file, ic_file=args.ic_file_name)