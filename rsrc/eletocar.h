void eletocar(double ell[6], double GM, double pos[3], double vit[3])
{
  double a,k,h,q,p,lambda,n;
  double fle,asr,lf,phi,psi,tp,tq,dg;
  double x1,y1,vx1,vy1,corf,sam1,cis2;
//  cout << "entree dans eletocar " << endl;
  a=ell[0];
  k=ell[1];
  h=ell[2];
  q=ell[3];
  p=ell[4];
  lambda=ell[5];
  n=sqrt(GM/(a*a*a));
  fle=lambda-k*sin(lambda)+h*cos(lambda);
  corf=(lambda-fle+k*sin(fle)-h*cos(fle))/(1.-k*cos(fle)-h*sin(fle));
  double shift=1.e-8;
  while (fabs(corf)>shift) {
    fle+=corf;
    corf=(lambda-fle+k*sin(fle)-h*cos(fle))/(1.-k*cos(fle)-h*sin(fle));
    shift *= 1.1;
//    cout << "corf " << corf << endl;
  }
  lf=-k*sin(fle)+h*cos(fle);
  sam1=-k*cos(fle)-h*sin(fle);
  asr=1./(1.+sam1);
  phi=sqrt(1.-k*k-h*h);
  psi=1./(1.+phi);
  x1=a*(cos(fle)-k-psi*h*lf);
  y1=a*(sin(fle)-h+psi*k*lf);
  vx1=n*asr*a*(-sin(fle)-psi*h*sam1);
  vy1=n*asr*a*(cos(fle)+psi*k*sam1);
  cis2=2.*sqrt(1.-p*p-q*q);
  tp=1.-2.*p*p;
  tq=1.-2.*q*q;
  dg=2.*p*q;
  pos[0]=x1*tp+y1*dg;
  pos[1]=x1*dg+y1*tq;
  pos[2]=(-x1*p+y1*q)*cis2;
  vit[0]=vx1*tp+vy1*dg;
  vit[1]=vx1*dg+vy1*tq;
  vit[2]=(-vx1*p+vy1*q)*cis2;
//  cout << "sortie de eletocar " << endl;
}

void eletocar(long double ell[6], long double GM, long double pos[3], long double vit[3])
{
  long double a,k,h,q,p,lambda,n;
  long double fle,asr,lf,phi,psi,tp,tq,dg;
  long double x1,y1,vx1,vy1,corf,sam1,cis2;
  a=ell[0];
  k=ell[1];
  h=ell[2];
  q=ell[3];
  p=ell[4];
  lambda=ell[5];
  n=sqrtl(GM/(a*a*a));
  fle=lambda-k*sinl(lambda)+h*cosl(lambda);
  corf=(lambda-fle+k*sinl(fle)-h*cosl(fle))/(1.L-k*cosl(fle)-h*sinl(fle));
  while (fabsl(corf)>1.e-14) {
    fle+=corf;
    corf=(lambda-fle+k*sinl(fle)-h*cosl(fle))/(1.L-k*cosl(fle)-h*sinl(fle));
//    cout << "corf " << corf << endl;
  }
  lf=-k*sinl(fle)+h*cosl(fle);
  sam1=-k*cosl(fle)-h*sinl(fle);
  asr=1.L/(1.L+sam1);
  phi=sqrtl(1.L-k*k-h*h);
  psi=1.L/(1.L+phi);
  x1=a*(cosl(fle)-k-psi*h*lf);
  y1=a*(sinl(fle)-h+psi*k*lf);
  vx1=n*asr*a*(-sinl(fle)-psi*h*sam1);
  vy1=n*asr*a*(cosl(fle)+psi*k*sam1);
  cis2=2.L*sqrtl(1.L-p*p-q*q);
  tp=1.L-2.L*p*p;
  tq=1.L-2.L*q*q;
  dg=2.L*p*q;
  pos[0]=x1*tp+y1*dg;
  pos[1]=x1*dg+y1*tq;
  pos[2]=(-x1*p+y1*q)*cis2;
  vit[0]=vx1*tp+vy1*dg;
  vit[1]=vx1*dg+vy1*tq;
  vit[2]=(-vx1*p+vy1*q)*cis2;
}
