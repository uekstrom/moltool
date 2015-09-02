#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <exception>
#include "taylor.h"
#include "taylor_newton.h"
#include "quaternion.h"

using namespace std;

int verbosity = 0; // 0 - standard amount of output, -1 no output, >0 debugging output.

bool unit_masses = false; //If true, give all element the same mass (1 au)

const char *element_syms[] = 
{ "X", "H", "He", 
  "Li", "Be", "B", "C", "N", "O", "F", "Ne",
  "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
  "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
  "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",

  "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
  "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",

  "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
  "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo",0};

//From http://en.wikipedia.org/wiki/List_of_elements_by_atomic_weight
const double atomic_weights[] = 
{ 0.0, 1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, 14.0067, 15.9994,
18.9984032, 20.1797, 22.98976928, 24.3050, 26.9815386, 28.0855, 30.973762,
32.065, 35.453, 39.948, 39.0983, 40.078, 44.955912, 47.867, 50.9415, 51.9961,
54.938045, 55.845, 58.6934, 58.933195, 63.546, 65.38, 69.723, 72.64,
74.92160, 78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585, 91.224, 92.90638,
95.96, 98, 101.07, 102.90550, 106.42, 107.8682, 112.411, 114.818, 
118.710, 121.760, 127.60, 126.90447, 131.293, 132.9054519, 137.327, 138.90547, 140.116,
140.90765, 144.242, 145, 150.36, 151.964, 157.25, 158.92535, 162.500, 164.93032, 167.259, 
168.93421, 173.054, 174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 
196.966569, 200.59, 204.3833, 207.2, 208.98040, 210, 210, 222, 223, 226, 227, 232.03806, 
231.03588, 238.02891, 237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262, 261, 262, 
  266, 264, 267, 268, 271, 272, 285, 284, 289, 288, 292, 294 /* element 118 */ };

int atomic_number(const string &sym)
{
    for (int i=0;element_syms[i];i++)
	if (sym == element_syms[i])
	    return i;
    return -1;
}

double atomic_mass_amu(int Z)
{
    assert(Z>=0);
    assert(Z<=118);
    if (unit_masses)
	return 1;
    else
	return atomic_weights[Z];
}


void read_xyz(vector<string> &syms, vector<vec<double,3> > &r, const string &path)
{
    istream *psrc;
    if (path == "-")
	psrc = &cin;
    else
	psrc = new ifstream(path.c_str());
    istream &src = *psrc;
    //Skip nr atoms and comment line
    int natoms;
    src >> natoms;
    if (verbosity > 0)
	cerr << "Reading " << natoms << " atoms from " << path << endl;
    src.ignore(10000,'\n');
    src.ignore(10000,'\n');
    for (int i=0;i<natoms;i++)
    {
	string sym;
	vec<double,3> pos;
	src >> sym;
	src >> pos[0] >> pos[1] >> pos[2];
	syms.push_back(sym);
	r.push_back(pos);
    }
    if (psrc != &cin)
	delete psrc;
}

void write_xyz(const vector<string> &sym, const vector<vec<double,3> > &r, ostream &dst)
{
    dst.precision(15);
    dst << fixed;
    dst << r.size() << endl;
    dst << endl;
    for (int i = 0;i<r.size();i++)
	dst << sym[i] << "    " << r[i][0] << " " << r[i][1] << " " << r[i][2] << endl;    
}


vec<double,3> center_of_mass(const vector<string> &sym, const vector<vec<double,3> > &r)
{
    vec<double,3> r_cm = 0;
    double mass = 0;
    for (int i = 0;i<sym.size();i++)
    {
	int Z = atomic_number(sym[i]);
	double m;
	if (Z<0)
	{
	    if (verbosity>=0)
		cerr << "Unknown element " << sym[i] << 
		    ", assuming mass = 1 amu" << endl;
	    m = 1;
	}
	else
	{
	    m = atomic_mass_amu(Z);
	}
	r_cm += m*r[i];
	mass += m;
    }
    return (1/mass)*r_cm;
}

template<class T>
T rotated_square_dev(const vector<vec<double,3> > &q1, 
		     const vector<vec<double,3> > &q2,
		     const quaternion<T> &R)
{
    T dev = 0;
    vec<T,3> r1,r2,rr;
    assert(q1.size() == q2.size());
    for (int i=0;i<q1.size();i++)
    {
	r1 = q1[i];
	r2 = q2[i];
	R.rotate(rr,r2);	    
	dev += (r1-rr).abs2();
    }
    return dev;
}

void rotate_geom(vector< vec<double,3> > &pos, const quaternion<double> &R)
{
    for (int i = 0;i<pos.size();i++)
    {
	vec<double,3> tmp;
	R.rotate(tmp,pos[i]);
	pos[i] = tmp;
    }
}

// Return a second order expansion of exp(i*x + j*y + k*z) at
// x=y=z=0
quaternion< taylor<3,2> > rotor_expansion(void)
{
    quaternion< taylor<3,2> > q;
    q[0] = 1;
    q[0][4] = -0.5;
    q[0][7] = -0.5;
    q[0][9] = -0.5;
    for (int i=0;i<3;i++)
	q[i+1] = taylor<3,2>(0,i);
    return q;
}


// Read a molecule from path argv[0], and an axis from argv[1]..argv[3]
int rotate(int argc, char *argv[], ostream &dst)
{
    if (argc < 4)
    {
	cerr << "Usage: moltool -r MOLECULE RX RY RZ [ANGLE]\n";
	return -1;
    }
    vector<string> sym;
    vector<vec<double,3> > r;
    read_xyz(sym,r,argv[0]);
    vec<double,3> axis;
    for (int i=0;i<3;i++)
    {
	stringstream ss;
	ss << argv[i+1];
	ss >> axis[i];
	if (ss.fail())
	{
	    cerr << "Error reading rotation axis " << i << endl;
	    return -1;
	}
    }
    if (argc == 5)
    {
	double angle;
	stringstream ss;
	ss << argv[4];
	ss >> angle;
	// Axis-angle format
	axis = M_PI*angle/360*axis.normalized();
    }
    else
    {
	axis *= 0.5;
    }
    quaternion<double> R = quaternion_rotor(axis);
    rotate_geom(r,R);
    write_xyz(sym,r,dst);
    return 0;
}

int translate(int argc, char *argv[], ostream &dst)
{
    if (argc < 4)
    {
	cerr << "Usage: moltool -t MOLECULE X Y Z\n";
	return -1;
    }
    vector<string> sym;
    vector<vec<double,3> > r;
    read_xyz(sym,r,argv[0]);
    vec<double,3> d;
    for (int i=0;i<3;i++)
    {
	stringstream ss;
	ss << argv[i+1];
	ss >> d[i];
	if (ss.fail())
	{
	    cerr << "Error reading translation vector component " << i << endl;
	    return -1;
	}
    }
    for (int i=0;i<r.size();i++)
	r[i] += d;
    write_xyz(sym,r,dst);
    return 0;
}

int find_center_of_mass(int argc, char *argv[], ostream &dst)
{
    if (argc < 1)
    {
	cerr << "Usage: moltool -cm MOLECULE\n";
	return -1;
    }
    vector<string> sym;
    vector<vec<double,3> > r;
    read_xyz(sym,r,argv[0]);
    vec<double,3> cm = center_of_mass(sym,r);
    dst << cm[0] << " " << cm[1] << " " << cm[2] << endl;
    return 0;
}

int rot_fit_rms(int argc, char *argv[], ostream &dst)
{
    if (argc < 2)
    {
	cerr << 
	    "Usage: moltool -rms MOLECULE1 MOLECULE2\n"
	    " Translate MOLECULE2 to have the same center of mass\n"
	    " as MOLECULE1, then rotate MOLECULE2 around its c.m. to\n"
	    " match MOLECULE1 in the best possible (RMS) way.\n";
	return -1;
    }
    vector<string> refsym, rotsym;
    vector<vec<double,3> > refr,rotr;
    read_xyz(refsym,refr,argv[0]);
    read_xyz(rotsym,rotr,argv[1]);
    if (refr.size() != rotr.size())
    {
	cerr << "ERROR in -rms, molecule must have the same number of atoms.\n";
	return -1;
    }
    // Move both molecules to have their center of mass at the origin
    vec<double,3> ref_cm = center_of_mass(refsym,refr), 
	          rot_cm = center_of_mass(rotsym,rotr);
    if (verbosity > 0)
      {
	cerr << "Molecule 1 center of mass: " 
	     << ref_cm[0] << " " << ref_cm[1] << " " << ref_cm[2] << endl;
	cerr << "Molecule 2 center of mass: " 
	     << rot_cm[0] << " " << rot_cm[1] << " " << rot_cm[2] << endl;
      }  
    for (int i=0;i<refr.size();i++)
    {
	refr[i] -= ref_cm;
	rotr[i] -= rot_cm;
    }

    quaternion<double> Rtot = 1;
    double best = 1e100;
    double max_step = 1;
    for (int ii=0;ii<1000;ii++)
    {
	quaternion< taylor<3,2> > Rt;
	taylor<3,2> devt;
	Rt = rotor_expansion();
	devt = rotated_square_dev(refr,rotr,Rt);
	if (devt[0] > best + 1e-15)
	{
	    cerr << "Target function goes up, quitting.\n";
	    return -1;
	}
	best = devt[0];
	vec<double,3> x;
	if (levenberg_step(x,devt,max_step) != 0)
	{
	    cerr << "Micro-optimization not converging, quitting" << endl;
	    return -1;
	}
	double step_norm = sqrt(x.abs2());
	if (verbosity > 0)
	    cerr << "RMS: " << sqrt(best/rotr.size()) << 
		", step norm: " << step_norm << endl;
	quaternion<double> R = quaternion_rotor(x);
	rotate_geom(rotr,R);
	Rtot = R*Rtot;
	if (step_norm < 1e-8)
	    break;
    }
    
    if (verbosity>=0)
    {
	Rtot = log(Rtot);
	Rtot *=2;
	cerr.precision(15);
	cerr << "Translate " << -rot_cm[0] << " " << -rot_cm[1] << " " << -rot_cm[2] << endl;
	cerr << "Rotation axis: " << Rtot[1] << " " << Rtot[2] << " " << Rtot[3] << endl;
	cerr << "Translate " << ref_cm[0] << " " << ref_cm[1] << " " << ref_cm[2] << endl;
	cerr << "Repeat rotation by:" << endl;
	cerr << "moltool -t - " 
	     << -rot_cm[0] << " " << -rot_cm[1] << " " << -rot_cm[2] 
	     << " | moltool -r - " 
	     << Rtot[1] << " " << Rtot[2] << " " << Rtot[3]
	     << " | moltool -t - "
	     << ref_cm[0] << " " << ref_cm[1] << " " << ref_cm[2] << endl;
    }
    // Move MOLECULE2 to the original position of MOLECULE1
    for (int i=0;i<refr.size();i++)
	rotr[i] += ref_cm;
    write_xyz(rotsym,rotr,dst);
    return 0;
}

int main(int argc, char *argv[])
{
    for (int i=1;i<argc;i++)
    {
	if (argv[i] == string("-r") or argv[i] == string("--rotate"))
	    return rotate(argc-i-1,argv+i+1,cout);
	if (argv[i] == string("-t") or argv[i] == string("--translate"))
	    return translate(argc-i-1,argv+i+1,cout);
	if (argv[i] == string("-cm") or argv[i] == string("--center-of-mass"))
	    return find_center_of_mass(argc-i-1,argv+i+1,cout);
	else if (argv[i] == string("-rms"))
	    return rot_fit_rms(argc-i-1,argv+i+1,cout);
	else if (argv[i] == string("-v"))
	    verbosity++;
	else if (argv[i] == string("-q"))
	    verbosity--;
	else if (argv[i] == string("--unit-masses"))
	    unit_masses = true;
	else
	{
	    cerr << "Unknown option '" << argv[i] << "', quitting." << endl;
	    return -1;
	}
    }
    return 0;
}
