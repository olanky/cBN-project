#include <iostream>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <cmath>

// This code is to bombard atoms or molecule toward the target
using namespace std;

// Maximum number of possible atoms
const int N = 10000;	

int main(int argc, const char* argv[])
{
    string line;
    ifstream infile;
    ofstream ofile1,ofile2;
    
    char element[N],str[50];
    // Particle type added in new bombardment
    // NH: 1; N: 2; H: 3
    int i,j,n[N],newSpeciesType,num=0,num_final;
    double x[N],y[N],z[N],vx[N],vy[N],vz[N], cmx=0, cmy=0, cmz=0;
    // Pre-factor for the calculation of particle velocities:
    // vel_particle_factor = sqrt(3*k/mass) 
    // NH: 4.0757e-1; N: 4.2198e-1; H: 1.5730 
    double vel_particle,vel_particle_factor[3]={4.0757e-1,4.2198e-1,1.5730},theta,phi,theta2,phi2;

    strcpy (str,"MD iter: ");
    strcat (str,argv[2]);
    strcat (str,"000");
   
  //////////////
  // Usage hint
    if (argc < 5) {
        cerr<<"Usage:"<<argv[0]<<" [ID of bombardment]"<<" [time per bombardment (ps)]"<<" [temperature (K)]" << " [type of particles]" <<endl;
	cerr<<"Type: 1: NH; 2: N; 3: H" << endl;
        cerr<<"It will generate input file named bnnt.gen and veloc.dat"<<endl;
        return 1;
    }

    infile.open("./geo_stp.xyz");
    ofile1.open("./bnnt.gen");
    ofile2.open("./veloc.dat");
  
  //////////////  
  // Check the last frame of last run
    do {
        getline(infile,line);
    }while(line!=str);
  // Read in the coordinates and velocities
    while(infile){
        getline(infile, line);
        stringstream ss(line);
        ss >> element[num] >> x[num] >> y[num] >> z[num] >> vx[num] >> vy[num] >> vz[num];
        ++num;
    };
    --num; // This is a correction of the total number of atoms. (some side effect of while loop)
  
  //////////////////
  // Calculate the center of mass
    int num_B=0,num_N=0,num_H=0;
    num_final=num;
    for (i=0;i<num;i++){
      // Check the isolated atoms; make a mark on the those atoms and correct the actual atom number
        if ((x[i]*x[i]+y[i]*y[i]+z[i]*z[i])>1300) {
            element[i]='X';
            --num_final;
        }
        switch (element[i]){
            case 'B':
                cmx+=10.811*x[i];
                cmy+=10.811*y[i];
                cmz+=10.811*z[i];
                num_B++;
                break;
            case 'N':
                cmx+=14.007*x[i];
                cmy+=14.007*y[i];
                cmz+=14.007*z[i];
                num_N++;
                break;
	    case 'H':
		cmx+=1.008*x[i];
		cmy+=1.008*y[i];
		cmz+=1.008*z[i];
		num_H++;
		break;
            default:
                break;
        }
    }
    cmx=cmx/(num_B*10.811+num_N*14.007+num_H*1.008);
    cmy=cmy/(num_B*10.811+num_N*14.007+num_H*1.008);
    cmz=cmz/(num_B*10.811+num_N*14.007+num_H*1.008);
  
  //////////////////
  // Make a shift of the stucture to make sure cm is at (0,0,0)
    for (i=0,j=0;i<num;i++){
        x[i]-=cmx;
        y[i]-=cmy;
        z[i]-=cmz;
      // Double check the isolated atoms beyond the spherical boundary of 3 nm radius from center of mass
        if ((x[i]*x[i]+y[i]*y[i]+z[i]*z[i])>900) {
            element[i]='X';
            --num_final;
        }
        if(element[i]!='X'){
		n[i]=j;
		++j;
	}
    }
   
    newSpeciesType=atoi(argv[4]); 
  ///////////////
  // Manually generate random number for particle bombardment with even distribution 
    theta=(atoi(argv[1])%11+1)*M_PI/12;
    phi=(atoi(argv[1])%5)*2*M_PI/5;
    theta2=(atoi(argv[1])%7+1)*M_PI/8;
    phi2=(atoi(argv[1])%13)*2*M_PI/13;
  // The velocity corresponds to the set temperature
    vel_particle=vel_particle_factor[newSpeciesType-1]*sqrt(atoi(argv[3])); 
    
  ///////////////
  // format the output file
    switch (newSpeciesType) {
	case 1:
            ofile1 << num_final + 2 << "  C" << endl;
            ofile1 << "B  N  H" << endl;
	    break;
	case 2:
            ofile1 << num_final + 1 << "  C" << endl;
            ofile1 << "B  N" << endl;
	    break;
	case 3:
            ofile1 << num_final + 1 << "  C" << endl;
            ofile1 << "B  N  H" << endl;
	    break;
        default:
            break;
    }
    for (i=0; i<num; i++) {
        switch (element[i]) {
            case 'B':
                ofile1 << setw(4) << n[i]+1 << setw(4)<< "1" << "\t" << setw(12) << x[i] << "\t" << setw(12) << y[i] << "\t" << setw(12) << z[i] << endl;
                ofile2 << setw(12) << vx[i] << "\t" << setw(12) << vy[i] << "\t" << setw(12) << vz[i] <<endl;
                break;
            case 'N':
                ofile1 << setw(4) << n[i]+1 << setw(4)<< "2" << "\t" << setw(12) << x[i] << "\t" << setw(12) << y[i] << "\t" << setw(12) << z[i] << endl;
                ofile2 << setw(12) << vx[i] << "\t" << setw(12) << vy[i] << "\t" << setw(12) << vz[i] <<endl;
                break;
	    case 'H':
                ofile1 << setw(4) << n[i]+1 << setw(4)<< "3" << "\t" << setw(12) << x[i] << "\t" << setw(12) << y[i] << "\t" << setw(12) << z[i] << endl;
                ofile2 << setw(12) << vx[i] << "\t" << setw(12) << vy[i] << "\t" << setw(12) << vz[i] <<endl;
                break;
            default:
                break;
        }
    }
    switch (newSpeciesType) {
	case 1:
    ofile1 << setw(4) << num_final+1 << setw(4)<< "2" << "\t" << setw(12) << 30*sin(theta)*cos(phi) << "\t" << setw(12) << 30*sin(theta)*sin(phi) << "\t" << setw(12) << 30*cos(theta) << endl;
    ofile1 << setw(4) << num_final+2 << setw(4)<< "3" << "\t" << setw(12) << 30*sin(theta)*cos(phi)+1.07*sin(theta2)*cos(phi2) << "\t" << setw(12) << 30*sin(theta)*sin(phi)+1.07*sin(theta2)*sin(phi2) << "\t" << setw(12) << 30*cos(theta)+1.07*cos(theta2) << endl;
    ofile2 << setw(12) << -vel_particle*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_particle*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_particle*cos(theta) <<endl;
    ofile2 << setw(12) << -vel_particle*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_particle*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_particle*cos(theta) <<endl;
	    break;
	case 2:
	    ofile1 << setw(4) << num_final+1 << setw(4)<< "2" << "\t" << setw(12) << 30*sin(theta)*cos(phi) << "\t" << setw(12) << 30*sin(theta)*sin(phi) << "\t" << setw(12) << 30*cos(theta) << endl;
	    ofile2 << setw(12) << -vel_particle*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_particle*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_particle*cos(theta) <<endl;
	    break;
	case 3:
            ofile1 << setw(4) << num_final+1 << setw(4)<< "3" << "\t" << setw(12) << 30*sin(theta)*cos(phi) << "\t" << setw(12) << 30*sin(theta)*sin(phi) << "\t" << setw(12) << 30*cos(theta) << endl;
            ofile2 << setw(12) << -vel_particle*sin(theta)*cos(phi) << "\t" << setw(12) << -vel_particle*sin(theta)*sin(phi) << "\t" << setw(12) << -vel_particle*cos(theta) <<endl;
	    break;
        default:
            break;
    }
    
    infile.close();
    ofile1.close();
    ofile2.close();
    return 0;
}
