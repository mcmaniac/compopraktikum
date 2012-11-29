// Elektronenstreuung in Materie
// Zu ersetzende oder fehlende Programmteile sind gekennzeichnet

// Durch Kompilieren des Programms kann die Ausführungszeit gewaltig verkürzt werden (s. ROOT-Manual)

#include<iostream>
#include<fstream>
#include<math.h>

// constants

const double pi = 3.14159265358;      // pi
const double me = 0.511;              // electron mass [MeV]

const int N_ELEC = 50;           // number of electrons
const double X_MAX = 100;        // x-size of the detector
const double Y_MAX = 100;        // y-size of the detector
const double E_0 = 20;           // energy of electrons
const double E_MIN = 0.02;       // minimum energy of electrons

class electron: public TObject {
public :
    electron(double x1, double y1, double z1, double E1 , double theta1, double phi1);
    ~electron();
    double x;
    double y;
    double z;
    double E;
    double theta;
    double phi;
};

electron::electron(double x1, double y1, double z1, double E1 , double theta1, double phi1) {
    x=x1; y=y1; z=z1; E=E1; theta=theta1; phi=phi1;
}

electron::~electron();

int nextcoord(electron *e1, electron *e2) 
{

// Keep starting position of struck electron e2

    e2->x = e1->x;
    e2->y = e1->y;
    e2->z = e1->z;
    double otheta = e1->theta;  // old theta
    double ophi   = e1->phi;    // old phi
    
// calculate new coordinates

// Moeller scattering
    double n_theta = f_theta->GetRandom();
    double n_phi = f_phi->GetRandom();
    double len= f_length->GetRandom()*e1->E/10; // scattering length (scales with E)
    double p = sqrt((e1->E+me)**2-me**2);      // Momentum of e1 in target system

    // momentum
    double p_a = p;
    double p_A = sqrt( pow(e1->E + me, 2) - pow(me, 2) );
    // invariant mass
    double s = pow(p_a + p_A, 2);
    // lorentz-factor
    double gamma = (e1->E + me) / sqrt(s);

// transform angles into target system (fixed)
    double theta1 = atan( sin(n_theta) / (gamma * (1 + cos(n_theta))) );
    double theta2 = atan( sin(n_theta) / (gamma * (1 - cos(n_theta))) );

// calculate new momenta and energies (fixed)
    double p1 = ( 2 * me * (e1->E + me) * p_a * cos(theta1) )
              / ( pow(e1->E + me, 2) - pow(p_a * cos(theta1), 2) );
    double p1 = ( 2 * me * (e1->E + me) * p_a * cos(theta2) )
              / ( pow(e1->E + me, 2) - pow(p_a * cos(theta2), 2) );

// stimmt!
    e1->E = (sqrt(p1**2+me**2)-me)-0.05*len/e1->E;
    e2->E = (sqrt(p2**2+me**2)-me)-0.05*len/e1->E;

// stop if no energy left
    if (e1->E <= E_MIN) {return 0;}

// transform into cartesian coordinates (x',y',z') 
    double x_ = len * sin(theta1) * cos(n_phi);
    double y_ = len * sin(theta1) * sin(n_phi);
    double z_ = len * cos(theta1);
    
// Rotation of e1 relative to  previous direction of e1 
    TVector3 v_(x_, y_, z_);
    v_.RotateX((-1.0)*otheta)
    v_.RotateZ((-1.0)*ophi + pi / 2.0);
    x_ = v_.X();
    y_ = v_.Y();
    z_ = v_.Z();
    
// New position and direction of e1 (fixed)

    e1->x += v_.X();
    e1->y += v_.Y();
    e1->z += v_.Z();
    e1->theta = theta1;
    e1->phi = n_phi;

// Transformation of e2 relative to  previous direction of e1 (fixed)
    double len_ = f_length->GetRandom() * e2->E / 10.0;
    double x__ = len_ * sin(theta2) * cos(n_phi);
    double y__ = len_ * sin(theta2) * sin(n_phi);
    double z__ = len_ * cos(theta2);

    TVector3 v__(x__, y__, z__);
    v__.RotateX((-1.0)*otheta);
    v__.RotateZ((-1.0)*ophi + pi / 2.0);

// Direction of e2 (fixed)

    e2->x     = e2->x + x__;
    e2->y     = e2->y + y__;
    e2->z     = e2->z + z__;
    e2->theta = acos( e2->z / sqrt(pow(e2->x, 2) + pow(e2->y, 2) + pow(e2->z, 2)) );
    e2->phi   = atan2(e2->y, e2->x);


// stop if electron leaves detector
    if ((e1->x < 0) || (e1->x>X_MAX)) {cout << "detector left" << endl; return 0;} 
    if ((e1->y < 0) || (e1->y>Y_MAX)) {cout << "detector left" << endl; return 0;}

  // return value 1 means success
    return 1;
}

void electron_template(void)
{
    ofstream outfile;
    ofstream logfile;
    TList electrons;    // List of electron objects (see root User's guide)
    electron* e1;       // first electron (whose trajectory is followed)
    electron* e2;       // second electron (which is just produced and traced afterwards)
    
// define graph
    TView* view;
    TPolyLine3D* traj = NULL;
    TCanvas* c1 = new TCanvas("SimCanvas", "SimCanvas", 1);
    int k = 0;
    
    // define functions for random number generator (kann so bleiben)
    TF1 f_length("f_length", "exp(-x)" ,0 ,100);   // mean free path
    TF1 f_phi("f_phi", "1" ,0 , 2*pi);             // azimuthal angle
    TF1 f_theta("f_theta",                         // Moeller scattering in c.m. system
		"(1+cos(x/2)**4)/sin(x/2)**4+2/(sin(x/2)*cos(x/2))**2+(1+sin(x/2)**4)/cos(x/2)**4"
		,0.2 ,pi/2); // The cross section cannot be integrated to 0 degrees, 
    // a too small starting angle will never produce
    // branches, therefore 0.2 rad is a good compromise
    
    c1->cd();
    view = TView::CreateView(1);
    view->SetRange(0,0,0,100,100,100);
    view->ShowAxis();
    
    for (int n = 0; n < N_ELEC; n++) {
	// loop for electrons
	cout << "electron #" << n << endl;
	electrons.Add(new electron(50, 50 ,0 ,E_0, 0, 0));  // push new electron on stack
	while (1) {
	    e1 = (electron*)electrons.First();              // pop next starting point from stack
	    if (e1==NULL) break;
	    electrons.Remove(e1);                           // remove it from the stack  
	    traj = new TPolyLine3D(1024);
	    k = 0;
	    traj->SetPoint(k, e1->x, e1->y, e1->z);
	    traj->Draw("ogl");                                                                   
	    k++;
	    while (1) {
		e2 = new electron(0, 0 ,0 ,0 , 0, 0);
		int ret=nextcoord(e1, e2);
		traj->SetPoint(k, e1->x, e1->y, e1->z);
		traj->Draw("ogl");
		k++;
		if (ret == 0) break;
		if (e2->E > E_MIN) electrons.Add(e2);     // put secondary electron on the stack
	    }  // end while
      }
      c1->Modified();
      c1->Update();                                                                   
  }  // end for
}

