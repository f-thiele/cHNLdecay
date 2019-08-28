
#include "Config.h"
#include "HNL.h"
#include "Lepton.h"
#include "Logger.h"
#include "Meson.h"
#include "auxfunctions.h"
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>

#include <chrono>
#include <thread>



// Tables

Double_t get_fhhp_V(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.76;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.295;
	else if(hID == 531 and hpID == 433) f=0.95;
	else if(hID == 531 and hpID == 323) f=0.291;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_V" << std::endl;
	return f; 
}

Double_t get_fhhp_A0(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.69;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.231;
	else if(hID == 531 and hpID == 433) f=0.67;
	else if(hID == 531 and hpID == 323) f=0.289;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_A0" << std::endl;
	return f; 
}
Double_t get_fhhp_A1(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.66;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.269;
	else if(hID == 531 and hpID == 433) f=0.70;
	else if(hID == 531 and hpID == 323) f=0.287;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_A1" << std::endl;
	return f; 
}
Double_t get_fhhp_A2(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t f(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) f=0.62;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) f=0.282;
	else if(hID == 531 and hpID == 433) f=0.75;
	else if(hID == 531 and hpID == 323) f=0.286;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_A2" << std::endl;
	return f; 
}

Double_t get_sigmahhp_V(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=0.57;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=0.875;
	else if(hID == 531 and hpID == 433) sigma=0.372;
	else if(hID == 531 and hpID == 323) sigma=0.516; // changed
	else std::cerr << "ERROR: in get_sigmahhp_V" << std::endl;
	return sigma; 
}
Double_t get_sigmahhp_A0(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=0.59;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=0.796;
	else if(hID == 531 and hpID == 433) sigma=0.350;
	else if(hID == 531 and hpID == 323) sigma=-0.383;
	else std::cerr << "ERROR: in get_sigmahhp_A0" << std::endl;
	return sigma; 
}
Double_t get_sigmahhp_A1(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=0.78;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=0.54;
	else if(hID == 531 and hpID == 433) sigma=0.463;
	else if(hID == 531 and hpID == 323) sigma=0.;
	else std::cerr << "ERROR: in get_sigmahhp_A1" << std::endl;
	return sigma; 
}
Double_t get_sigmahhp_A2(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t sigma(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) sigma=1.40;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) sigma=1.34;
	else if(hID == 531 and hpID == 433) sigma=1.04;
	else if(hID == 531 and hpID == 323) sigma=1.05;
	else std::cerr << "ERROR: in get_sigmahhp_A2" << std::endl;
	return sigma; 
}
//
Double_t get_xihhp_V(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=0.;
	else if(hID == 531 and hpID == 433) xi=0.561;
	else if(hID == 531 and hpID == 323) xi=2.1;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}
Double_t get_xihhp_A0(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=0.055;
	else if(hID == 531 and hpID == 433) xi=0.600;
	else if(hID == 531 and hpID == 323) xi=1.58;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}
Double_t get_xihhp_A1(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=0.;
	else if(hID == 531 and hpID == 433) xi=0.510;
	else if(hID == 531 and hpID == 323) xi=1.06;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}
Double_t get_xihhp_A2(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t xi(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) xi=0.41;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) xi=-0.21;
	else if(hID == 531 and hpID == 433) xi=0.070;
	else if(hID == 531 and hpID == 323) xi=-0.074;
	else std::cerr << "ERROR: in get_xihhp_V" << std::endl;
	return xi; 
}


Double_t get_MhP(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t M(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) M=6275.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) M=5279.;
	else if(hID == 531 and hpID == 433) M=6275.;
	else if(hID == 531 and hpID == 323) M=5367.;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_MhP" << std::endl;
	return M; 
}
Double_t get_MhV(Meson h, Meson hp){
	int hID = h.getPdgId();
	int hpID = hp.getPdgId();
	Double_t M(0);
	
	if((hID == 521 and hpID == 423)or(hID == 511 and hpID == 413)) M=6331.;
	else if((hID == 521 and hpID == 113)or(hID == 511 and hpID == 213)) M=5325.;
	else if(hID == 531 and hpID == 433) M=6331.;
	else if(hID == 531 and hpID == 323) M=5415.;
	else std::cerr << "ERROR: hID = "<<hID<<", hpID = " << hpID << ", in get_fhhp_MhV" << std::endl;
	return M; 
}
 
 


// Functions

Double_t compute_FF_V(Meson h, Meson hp, Double_t q2){
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_V(h, hp);
	Double_t sigmahhp = get_fhhp_V(h, hp);
	Double_t xihhp = get_xihhp_V(h, hp);
	Double_t MP = get_MhP(h, hp);
	Double_t MV = get_MhV(h, hp);
	
	res = (fhhP)/((1-(q2/pow(MV,2)))*(1-sigmahhp*(q2/pow(MV,2)) - xihhp*pow(q2,2)/pow(MV,4)));
	//if((hID == 531 and hpID == 323))std::cout<<"V: "<<res<<std::endl;
	return res;
}

Double_t compute_FF_A0(Meson h, Meson hp, Double_t q2){
	
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_A0(h, hp);
	Double_t sigmahhp = get_fhhp_A0(h, hp);
	Double_t xihhp = get_xihhp_A0(h, hp);
	Double_t MP = get_MhP(h, hp);
	Double_t MV = get_MhV(h, hp);
	
	res = (fhhP)/((1-q2/pow(MP,2))*(1-sigmahhp*q2/pow(MV,2) - xihhp*pow(q2,2)/pow(MV,4)));
	//std::cout<<"A0: "<<res<<std::endl;
	return res;
}

Double_t compute_FF_A1(Meson h, Meson hp, Double_t q2){
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_A1(h, hp);
	Double_t sigmahhp = get_fhhp_A1(h, hp);
	Double_t xihhp = get_xihhp_A1(h, hp);
	Double_t MP = get_MhP(h, hp);
	Double_t MV = get_MhV(h, hp);
	
	res = (fhhP)/((1-sigmahhp*q2/pow(MV,2) - xihhp*pow(q2,2)/pow(MV,4)));
	//std::cout<<"A1: "<<res<<std::endl;
	return res;
}

Double_t compute_FF_A2(Meson h, Meson hp, Double_t q2){
	Double_t res;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	Double_t fhhP = get_fhhp_A2(h, hp);
	Double_t sigmahhp = get_fhhp_A2(h, hp);
	Double_t xihhp = get_xihhp_A2(h, hp);
	Double_t MP = get_MhP( h, hp);
	Double_t MV = get_MhV( h, hp);
	
	res = (fhhP)/((1-sigmahhp*q2/pow(MV,2) - xihhp*pow(q2,2)/pow(MV,4)));
	//std::cout<<"A2: "<<res<<std::endl;
	return res;
}


// -------------
Double_t compute_FF_g(Meson h, Meson hp, Double_t q2){
	Double_t res, V;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	V = compute_FF_V(h, hp, q2);
	res = V/(mh+mhp);
	return res;
}

Double_t compute_FF_aplus(Meson h, Meson hp, Double_t q2){
	Double_t res, A2;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	A2 = compute_FF_A2(h, hp, q2);
	res = -1*A2/(mh+mhp);
	return res;
}

Double_t compute_FF_f(Meson h, Meson hp, Double_t q2){
	Double_t res, A1;
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	A1 = compute_FF_A1(h, hp, q2);
	res = (mh+mhp)*A1;
	return res;
}

Double_t compute_FF_aminus(Meson h, Meson hp, Double_t q2){
	Double_t res, A0;
	Double_t aplus = compute_FF_aplus(h, hp, q2);
	Double_t f = compute_FF_f(h, hp, q2);
	size_t hID = h.getPdgId(); size_t mh = h.getMass(); 
	size_t hpID = hp.getPdgId(); size_t mhp = hp.getMass(); 
	A0 = compute_FF_A0(h, hp, q2);
	res = (2*mhp*A0-f-(pow(mh,2)-pow(mhp,2))*aplus)/q2;
	return res;
}



Double_t compute_ffactor(Meson h, Meson hp, Double_t q2, int opt) {
	
	// Uses equation (C.16) of arXiv paper
	
	Double_t 	Mpole, mh, mhp, a0, a1, a2, z, tplus, t0, factor, sum;
	Int_t 		IDh, IDhp;
	
	mh =  h.getMass(); 		mhp =  hp.getMass();
	IDh = h.getPdgId();		IDhp = hp.getPdgId();
	
	switch(opt){
		
		case 1:				//form factor fplus
				
			switch(IDhp){	
				case 411:			//D
				//case 413:			//Dstar
				case 421:			//D0bar
				//case 423:			//D0starbar
				case 431:			//Ds
				//case 433:			//Dsstar
					//Mpole = inf;
					a0 = 0.909;
					a1 = -7.11;
					a2 = 66;
					factor = 1.;
					break;
				case 321:			//K
				//case 323:			//Kstar
					Mpole = 5325;
					a0 = 0.360;
					a1 = -0.828;
					a2 = 1.1;
					factor = 1/(1-q2/pow(Mpole, 2));
					break;
				case 211:			//pi
				case 111:			//pi0
					Mpole = 5325;
					a0 = 0.404;
					a1 = -0.68;
					a2 = -0.86;	
					factor = 1/(1-q2/pow(Mpole, 2));
					break;
			}
			break;
		case 0: // form factor f0
			switch(IDhp){	
				case 411:			//D
				case 413:			//Dstar
				case 421:			//D0bar
				case 423:			//D0starbar
				case 431:			//Ds
				//  case 433:			//Dsstar
					//Mpole = inf;
					a0 = 0.794;
					a1 = -2.45;
					a2 = 33;
					factor = 1.;
					break;
				case 321:			//K
				//case 323:			//Kstar
					Mpole = 5650;
					a0 = 0.233;
					a1 = 0.197;
					a2 = 0.18;
					factor = 1/(1-q2/pow(Mpole, 2));
					//std::cout<<"factor " <<std::endl;
					break;
				case 211:			//pi
				case 111:			//pi0
					Mpole = 5650;
					a0 = 0.490;
					a1 = -1.61;
					a2 = 0.93;	
					factor = 1/(1-q2/pow(Mpole, 2));
					break;
					//std::cout<<"factor " <<std::endl;
			}
			break;
				
		case 2:		// g
			//std::cout<<"g"<<compute_FF_g(h, hp, q2)<<std::endl;
			return compute_FF_g(h, hp, q2);
			//return 1.;
			break;
		case 3:		// f	
			//std::cout<<"f:"<<compute_FF_f(h, hp, q2)<<std::endl;
			//return 1.;
			return compute_FF_f(h, hp, q2);
			break;
		case 4:		// aplus	
			//std::cout<<"aplus:"<<compute_FF_aplus(h, hp, q2)<<std::endl;
			//return 1.;
			return compute_FF_aplus(h, hp, q2);
			break;
		case 5:		// aminus	
			//return 1.;
			//std::cout<<"aminus:"<<compute_FF_aminus(h, hp, q2)<<std::endl;
			return compute_FF_aminus(h, hp, q2);
			break;

	}
	
	
	tplus = pow(mh+mhp, 2);
	t0 	  = (mh+mhp)*pow(sqrt(mh)-sqrt(mhp),2);
	
	z = sqrt(tplus-q2) - sqrt(tplus-t0);
	z = z/(sqrt(tplus-q2) + sqrt(tplus-t0));
	
	sum  = a0;
	sum += a1*(z 		- 1/3*pow(z,3));
	sum += a2*(pow(z,2)	+ 2/3*pow(z,3));
	
	
	return factor*sum;
}
