void temp(){

//px=0.297 py=-0.335 pz=0.968 -162.586 0.129
//
//31.241 0.008 -1247.014
//px=0.063 py=0.273 pz=0.417 -162.334 1.653 6497.973

double px = 0.063;
double py=0.273;
double pz=0.417;

double qx = 31.241;
double qy=0.008;
double qz=-1247.014;

double a = -tan(0.025);
double b=0;
double c=1.0;
double d=-1.0*tan(0.025)*6500.0*sin(0.025)-6500.0*cos(0.025);
double tDenom = a*(qx-px) + b*(qy-py) + c*(qz-pz);
if (tDenom == 0) cout<<"error"<<endl;
double t = - ( a*px + b*py + c*pz + d ) / tDenom;
cout<<tDenom<<endl;
cout<<px+t*(qx-px)<<endl;
cout<<py+t*(qy-py)<<endl;
cout<<pz+t*(qz-pz)<<endl;
}


