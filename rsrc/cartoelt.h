using namespace std;

void cartoelt(double pos[3], double vit[3], double GM, double ell[6])
{
  double rayon,v2,a,lambda,k,h,q,p;
  double gx,gy,gz,gg,tp,tq,dg,cis2;
  double psi,ach,ack,adg,det,sm1,sm2,cf,sf;
  double x1,y1,vx1,vy1,fle;
  rayon=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
  v2=vit[0]*vit[0]+vit[1]*vit[1]+vit[2]*vit[2];
  a=GM*rayon/(2.*GM-rayon*v2);
  gx=pos[1]*vit[2]-pos[2]*vit[1];
  gy=pos[2]*vit[0]-pos[0]*vit[2];
  gz=pos[0]*vit[1]-pos[1]*vit[0];
  gg=sqrt(gx*gx+gy*gy+gz*gz);
  cis2=sqrt(0.5*(1.+gz/gg));
  q=-gy/(2.*gg*cis2);
  p=gx/(2.*gg*cis2);
  tp=1.-2.*p*p;
  tq=1.-2.*q*q;
  dg=2.*p*q;
  x1=tp*pos[0]+dg*pos[1]-2.*p*cis2*pos[2];
  y1=dg*pos[0]+tq*pos[1]+2.*q*cis2*pos[2];
  vx1=tp*vit[0]+dg*vit[1]-2.*p*cis2*vit[2];
  vy1=dg*vit[0]+tq*vit[1]+2.*q*cis2*vit[2];
  k=gg*vy1/GM-x1/rayon;
  h=-gg*vx1/GM-y1/rayon;
  psi=1./(1.+sqrt(1.-k*k-h*h));
  ach=1.-psi*h*h;
  ack=1.-psi*k*k;
  adg=psi*h*k;
  det=ach*ack-adg*adg;
  sm1=x1/a+k;
  sm2=y1/a+h;
  cf=(sm1*ack-sm2*adg)/det;
  sf=(ach*sm2-adg*sm1)/det;
  fle=atan2(sf,cf);
  lambda=fle-k*sf+h*cf;
  ell[0]=a;
  ell[1]=k;
  ell[2]=h;
  ell[3]=q;
  ell[4]=p;
  ell[5]=lambda;
}


void cartoelt(long double pos[3], long double vit[3], long double GM, long double ell[6])
{
  long double rayon,v2,a,lambda,k,h,q,p;
  long double gx,gy,gz,gg,tp,tq,dg,cis2;
  long double psi,ach,ack,adg,det,sm1,sm2,cf,sf;
  long double x1,y1,vx1,vy1,fle;
  rayon=sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
  v2=vit[0]*vit[0]+vit[1]*vit[1]+vit[2]*vit[2];
  a=GM*rayon/(2.*GM-rayon*v2);
  gx=pos[1]*vit[2]-pos[2]*vit[1];
  gy=pos[2]*vit[0]-pos[0]*vit[2];
  gz=pos[0]*vit[1]-pos[1]*vit[0];
  gg=sqrt(gx*gx+gy*gy+gz*gz);
  cis2=sqrt(0.5*(1.+gz/gg));
  q=-gy/(2.*gg*cis2);
  p=gx/(2.*gg*cis2);
  tp=1.-2.*p*p;
  tq=1.-2.*q*q;
  dg=2.*p*q;
  x1=tp*pos[0]+dg*pos[1]-2.*p*cis2*pos[2];
  y1=dg*pos[0]+tq*pos[1]+2.*q*cis2*pos[2];
  vx1=tp*vit[0]+dg*vit[1]-2.*p*cis2*vit[2];
  vy1=dg*vit[0]+tq*vit[1]+2.*q*cis2*vit[2];
  k=gg*vy1/GM-x1/rayon;
  h=-gg*vx1/GM-y1/rayon;
  psi=1./(1.+sqrt(1.-k*k-h*h));
  ach=1.-psi*h*h;
  ack=1.-psi*k*k;
  adg=psi*h*k;
  det=ach*ack-adg*adg;
  sm1=x1/a+k;
  sm2=y1/a+h;
  cf=(sm1*ack-sm2*adg)/det;
  sf=(ach*sm2-adg*sm1)/det;
  fle=atan2(sf,cf);
  lambda=fle-k*sf+h*cf;
  ell[0]=a;
  ell[1]=k;
  ell[2]=h;
  ell[3]=q;
  ell[4]=p;
  ell[5]=lambda;
}
