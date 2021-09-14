
using namespace std;
float pi = TMath::Pi();
float cT = .0789;
float mlam = 1.115;
float xmpr = 0.93827;
float xmk = 0.493677;

//float Elam, Ek, pk;
//float cthK, cthL, keL, pL, cthcm, angK;

float path_output(float Egamma, float pMin, float pMax);

float Path_Length(float egam, TF1* xs, float lowerMomentum, float upperMomentum);

float CthCM_Kaon(float egam, float p);

float CthCM(float egam, float p);

float CthL(float egam, float p);

float EgamMin(Float_t P_lambda);

float Emax(Float_t P_lambda);

float Pmax(Float_t Egam);

float p0(float egam);

void new_pl_calc(){

	ofstream plFile;
	plFile.open("pathLength_new.txt");
	
	//const int n = 16;
	const int n = 32;
	
	int binNumber = 69;
	float Luminosity = 0;
	float min, max;
	
	float Egam_list[]= {1.194, 1.219, 1.244, 1.269, 1.294, 1.32, 1.345, 1.37, 1.395, 1.421, 1.446, 1.471, 1.496, 1.521, 1.546, 1.571, 1.596, 1.621, 1.647, 1.672, 1.698, 1.723, 1.748, 1.774, 1.799, 1.824, 1.849, 1.875, 1.9, 1.925, 1.95, 1.975, 2.001, 2.026, 2.051, 2.076, 2.101, 2.126, 2.151, 2.177, 2.202, 2.227, 2.253, 2.278, 2.303, 2.328, 2.353, 2.378};

	float x[n], y[n];
	for (int momBin = 0; momBin < n; momBin++){
		min = 0.6 + 0.1*momBin;
		max = 0.7 + 0.1*momBin;
		
		int count = 0;
		float lc;
		float pl =0;
		for(int eBin=0; eBin<=46; eBin++){		
			if (Egam_list[eBin] >= EgamMin(min+.05)){

				pl += path_output(Egam_list[eBin], min, max);
				count++;
			}	
		}
	
		plFile << pl/count << endl;
		
		x[momBin] = min+.05;
		y[momBin] = pl/count;
		
	}

	plFile.close();
	
	
	
	/**************************** Plotting Path Length  ****************************/ 
	TCanvas *c1 = new TCanvas("c1","c1");
	TGraph *gr1 = new TGraph (n,x,y);	
	gr1->SetTitle("");
	gr1->GetXaxis()->SetTitle("P_{#Lambda Lab}");
	gr1->GetYaxis()->SetTitle("Path Length");
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	gr1->GetXaxis()->SetTitleSize(0.05);
	gr1->GetYaxis()->SetTitleSize(0.05);
	gr1->GetYaxis()->SetTitleOffset(1.0);
	gr1->GetYaxis()->SetTitleOffset(0.8);
	gr1->GetXaxis()->SetTitleOffset(0.8);

	gr1->SetLineWidth(2);
	gr1->SetMarkerStyle(kFullSquare);
	gr1->SetMarkerColor(kBlue);
	gr1->SetLineColor(kBlue);

	gr1->Draw("AP");
}


float path_output(float Egamma, float pMin, float pMax){
	
	ifstream file ("KLambda_legendre.txt");
	float s = 2*xmpr*Egamma+xmpr*xmpr;
	float pcm = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	
	gSystem->Load("libMathMore");
	string line;
	float Index, E, W, Order, C0, C1, C2, C3, C4, sig0, sig1, sig2, sig3, sig4;
	
	
	TF1* L;
	int count_line = 0;
		
	for(int i =0;i<47;i++){
		file >> Index >> E >> W >> Order >> C0 >> C1 >> C2 >> C3 >> C4 >> sig0 >> sig1 >> sig2 >> sig3 >> sig4;
		
		//cout << "Energy: " << E << endl;
		//Defining cross section
		TF1* XS = new TF1("Cross Section", "[0]*(1+ [1]*x + [2]*(0.5*(3*std::pow(x,2) - 1) ) + [3]*(0.5*(5*std::pow(x,3) - 3*x) ) + [4]*(0.125*(35*std::pow(x,4)-30*x*x+3) ))",-1,1);
		XS->SetParameters(C0, C1, C2, C3, C4, 1, 2, 3, 4);
		//L[lineno - 69] = XS;
		
		if (Egamma == E){
			L = XS;
		}
	}

	
	return Path_Length(Egamma, L, pMin, pMax);
}


float Path_Length(float egam, TF1* xs, float lowerMomentum, float upperMomentum){
	int totalLambda = 1000;
	float step = 0.0001;
	float X, Y, Z;
	int lambda = 0;
	float xAverage = 0;
	float targetLength = .40;
	float targetWidth = .04;
	
	int count_high = 0;
	int count_low = 0;
	while (lambda < totalLambda){
	
		float p = (upperMomentum - lowerMomentum)*(double)rand()/(double)(RAND_MAX) + lowerMomentum;
		
		float xs_check = ((float)rand()/(float)(RAND_MAX)) * 3.0;
		
		//cout << "cos(theta) =  " << CthCM(egam, p) << endl;
		float cross_section = xs->Eval(CthCM(egam, p));
		//float cross_section = xs->Eval(CthCM_Kaon(egam, p));
		float theta = acos( CthL(egam, p) );
		float upperLimit = CthCM(egam, upperMomentum);
		float lowerLimit = CthCM(egam, lowerMomentum);

		if (xs_check <= cross_section){
		
			float distanceTraveled = 0;
	
			float start = ((float)rand()/(float)(RAND_MAX)) * targetLength;
	
			float r = ((float)rand()/(float)(RAND_MAX)) *.02 + (-.01);

			float phi = (((float)rand()/(float)(RAND_MAX)) *2*pi + (-pi));
			
			float phi_start = (((float)rand()/(float)(RAND_MAX)) *2*pi + (-pi));
			float xInit = r*cos(phi_start);
			float yInit = r*sin(phi_start);
			float bInit = sqrt(xInit*xInit+yInit*yInit);
			
			float X, Y, Z, dx, dy, dz;			
			float prob = 1.0;
			float survivalValue = 0.0;
			X = xInit;
			Y = yInit;
			Z = start;
			
			//cout << "X,Y: " << X*100 << "," << Y*100 << endl;

			dx = step*cos(phi)*sin(theta);
			dy = step*sin(phi)*sin(theta);
			dz = step*cos(theta);
			while (survivalValue <= prob){
			
				distanceTraveled = distanceTraveled + step;
				
				X = X+dx;
				Y = Y+dy;
				Z = Z+dz;
				
				prob = exp(-(mlam/p)*(step/cT));
				survivalValue = ((float)rand()/(float)(RAND_MAX)) * 1.0;
				
				if (sqrt(X*X+Y*Y) > (targetWidth/2) || Z > targetLength){
					break;
				}	
			
			}

			xAverage += sqrt((X-xInit)*(X-xInit) + (Y-yInit)*(Y-yInit) + (Z-start)*(Z-start));
			++lambda;
		} //xs check
	}

	float path_length = (xAverage/float(totalLambda))*100.0;  //m->cm
	return path_length;
}

float CthCM(float egam, float p){

	float const mlam=1.115683, mprot=0.938272, mkp=0.493677;

	int i,j;
	float ecm,ekpcm,elamcm,pkpcm,plamcm;
	float etotlab,elamlab,ekplab,tmandel;
	//float egamcm,costhcmKp[2];
	float egamcm,costhcmKp;

	// CM frame: equation from PDG sections 43.2 and 43.4.2

	ecm = sqrt( mprot*mprot + 2.0*egam*mprot );
	ekpcm = (ecm*ecm - mlam*mlam + mkp*mkp)/(2.0*ecm);
	elamcm = (ecm*ecm - mkp*mkp + mlam*mlam)/(2.0*ecm);
	//cout << ecm << " " << ekpcm << " " << elamcm << endl;
	pkpcm = sqrt( ekpcm*ekpcm - mkp*mkp );
	plamcm = sqrt( elamcm*elamcm - mlam*mlam );
	//cout << pkpcm << " " << plamcm << endl;

	// Lab frame
	etotlab = egam + mprot;
	
	elamlab = sqrt(p*p + mlam*mlam);
	ekplab = etotlab - elamlab;

	// Mandelstam variable from PDG section 43.5.1
	tmandel = mprot*mprot + mlam*mlam -2.0*mprot*elamlab;
	//cout << elamlab << " " << ekplab << " " << tmandel << endl;

	// convert Lab to CM: PDG equation (43.6)
	egamcm = egam*mprot/ecm;

	// tmandel is invariant: t_cm = t_lab
	costhcmKp = (tmandel - mkp*mkp + 2.0*egamcm*ekpcm)/(2.0*egamcm*pkpcm);
	
	return costhcmKp;

}

float CthL(float egam, float p){	
	float const mlam=1.115683, mprot=0.938272, mkp=0.493677;

	float ecm,ekpcm,elamcm,pkpcm,plamcm;
	float etotlab,elamlab,ekplab,tmandel;
	float egamcm,costhcmKp;
	float costhcmLam, costhlabLam;

	// CM frame: equation from PDG sections 43.2 and 43.4.2
	ecm = sqrt( mprot*mprot + 2.0*egam*mprot );
	ekpcm = (ecm*ecm - mlam*mlam + mkp*mkp)/(2.0*ecm);
	elamcm = (ecm*ecm - mkp*mkp + mlam*mlam)/(2.0*ecm);

	pkpcm = sqrt( ekpcm*ekpcm - mkp*mkp );
	plamcm = sqrt( elamcm*elamcm - mlam*mlam );

	// Lab frame
	etotlab = egam + mprot;

	elamlab = sqrt(p*p + mlam*mlam);
	ekplab = etotlab - elamlab;

	// Mandelstam variable from PDG section 43.5.1
	tmandel = mprot*mprot + mlam*mlam -2.0*mprot*elamlab;

	// convert Lab to CM: PDG equation (43.6)
	egamcm = egam*mprot/ecm;

	// tmandel is invariant: t_cm = t_lab
	costhcmKp = (tmandel - mkp*mkp + 2.0*egamcm*ekpcm)/(2.0*egamcm*pkpcm);
	
	costhcmLam = (tmandel - mlam*mlam + 2.0*egamcm*elamcm)/(2.0*egamcm*plamcm);	
	//costhlabLam = (tmandel + 2.0*mprot*elamlab - mlam*mlam - mprot*mprot)/(2.0*mprot*p);
	costhlabLam = cos( asin((plamcm/p)*sin(pi-acos(costhcmKp))) );

	if (costhlabLam <-1.0 || costhlabLam >1.0) { cout << "costhlabLam: " << costhlabLam << endl; }
	
	return costhlabLam;
}


float EgamMin(Float_t P_lambda){
	Float_t Emin = 0.890938 *(P_lambda) + 0.168662;
	//Float_t Emin = (sqrt(P_lambda*P_lambda+mlam*mlam)+xmk-xmpr);
	return Emin;
}

float Emax(Float_t P_lambda){
	float mlam = 1.115;
	float xmpr = 0.93827;
	float xmk = 0.493677;

	Float_t E = sqrt((P_lambda*P_lambda)+(mlam*mlam))+xmk-xmpr;
	
	return E;
}

float Pmax(Float_t Egam){
	float mlam = 1.115;
	float xmpr = 0.93827;
	float xmk = 0.493677;
	Float_t Elam = Egam +xmpr-xmk;
	Float_t Pmax = sqrt(Elam*Elam-mlam*mlam);
	
	return Pmax;
}


float p0(float egam){
	float pi = TMath::Pi();
	float mlam = 1.115;
	float xmpr = 0.93827;
	float xmk = 0.493677;
	
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt(((s-xmpr*xmpr)*(s-xmpr*xmpr))/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	float xi = log( (pcm+sqrt(xmpr*xmpr+pcm*pcm)/xmpr) );
	
	float Elam = sqrt(pcm_prime*pcm_prime +mlam*mlam)*cosh(xi);
	float p = sqrt(Elam*Elam-mlam*mlam);

	return p;
}

