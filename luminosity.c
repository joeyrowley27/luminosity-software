
using namespace std;
float pi = TMath::Pi();
float cT = .0789;
float mlam = 1.115;
float xmpr = 0.93827;
float xmk = 0.493677;

//float Elam, Ek, pk;
//float cthK, cthL, keL, pL, cthcm, angK;

float path_output(float Egamma, float pMin, float pMax);
float nlam_calc(float Egamma, float pMin, float pMax);
float Path_Length(float egam, TF1* xs, float lowerMomentum, float upperMomentum);
float CthCM_Kaon(float egam, float p);
float CthCM(float egam, float p);
float CthL(float egam, float p);
float EgamMin(Float_t P_lambda);
float p0(float egam);

void new_lum_calc(){

	ofstream lumFile;
	lumFile.open("luminosity_updated.txt");
	
	ofstream updatedFlux ("updatedFlux.txt");
	
	//ifstream fluxFile ("gflux_file.txt");
	ifstream fluxFile ("../outFile.txt");
	
	TFile *f1 = new TFile ("/home/joey/analysis/luminosity/Updated/NickCGenerator/GenerateLam/input/fluxhisto.root");
	TH1F *h1;
	h1 = (TH1F*)f1->Get("Hflux");
	
	const int n = 16;
	//const int n = 32;
	
	int binNumber = 69;
	float Luminosity = 0;
	float min, max;
	
	float Egam_list[]= {1.194, 1.219, 1.244, 1.269, 1.294, 1.32, 1.345, 1.37, 1.395, 1.421, 1.446, 1.471, 1.496, 1.521, 1.546, 1.571, 1.596, 1.621, 1.647, 1.672, 1.698, 1.723, 1.748, 1.774, 1.799, 1.824, 1.849, 1.875, 1.9, 1.925, 1.95, 1.975, 2.001, 2.026, 2.051, 2.076, 2.101, 2.126, 2.151, 2.177, 2.202, 2.227, 2.253, 2.278, 2.303, 2.328, 2.353, 2.378};
	
	float targetDensity = 6E23*0.071;

	float x[n], y[n];
	for (int momBin = 0; momBin < n; momBin++){
		min = 0.6 + 0.1*momBin;
		max = 0.7 + 0.1*momBin;
		//min = 0.60 + 0.05*momBin;
		//max = 0.65 + 0.05*momBin;
		
		int count = 0;
		float lc;
		float pl =0;
		float nlam =0;
		//for(int eBin=0; eBin<binNumber; eBin++){
		for(int eBin=0; eBin<=46; eBin++){
		//for(int eBin=startingE; eBin<=46; eBin++){

			/*
			if (TMath::IsNaN(nlam_calc(Egam_list[eBin], min, max))){
				count++;
			}else{
				pl += path_output(Egam_list[eBin], min, max);
				nlam += nlam_calc(Egam_list[eBin], min, max);
				count++;
			}
			*/
			
			if (Egam_list[eBin] > EgamMin(min+.05)){
			
				pl += path_output(Egam_list[eBin], min, max);
				nlam += nlam_calc(Egam_list[eBin], min, max);
				count++;
			}

		}
		
		float gflux = h1->Integral(h1->FindFixBin(EgamMin(min+.05)), h1->FindFixBin(2.4), "");
		updatedFlux << gflux << endl;
		
		float photonLum = gflux*(6E23*0.071)*40.0;

		lumFile << (nlam*targetDensity*photonLum)/count << endl;
		
		//cout << (nlam*targetDensity*photonLum)/count << endl;
		
		x[momBin] = min+.05;
		//x[momBin] = min+.005;
		y[momBin] = pl/count;
		
	}

	lumFile.close();
	fluxFile.close();
	
	
	/**************************** Plotting Path Length  ****************************/ 
	TCanvas *c1 = new TCanvas("c1","c1");
	//TGraphErrors *gr1 = new TGraphErrors(n,x,y, 0, error);
	TGraph *gr1 = new TGraph (n,x,y);
	//TGraphErrors *gr1 = new TGraphErrors (n,x,y, 0, error);
	
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
	
	
	/**************************** Plotting cos(Theta)_CM  ****************************/
	//TCanvas *c2 = new TCanvas("c1","c1");
	//Theta->Draw();
}


float path_output(float Egamma, float pMin, float pMax){
	
	//ifstream file ("KLambda_legendre.txt");
	ifstream file ("test_file.txt");
	float s = 2*xmpr*Egamma+xmpr*xmpr;
	float pcm = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	
	gSystem->Load("libMathMore");
	string line;
	float Index, E, W, Order, C0, C1, C2, C3, C4, sig0, sig1, sig2, sig3, sig4;
	
	
	TF1* L;
	for(int i =0;i<48;i++){
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

float nlam_calc(float Egamma, float pMin, float pMax){
	
	//ifstream file ("KLambda_legendre.txt");
	ifstream file ("test_file.txt");
	
	gSystem->Load("libMathMore");
	string line;
	float Index, E, W, Order, C0, C1, C2, C3, C4, sig0, sig1, sig2, sig3, sig4;
	
	
	TF1* L;
	int count_line = 0;
	for(int i =0;i<48;i++){
		file >> Index >> E >> W >> Order >> C0 >> C1 >> C2 >> C3 >> C4 >> sig0 >> sig1 >> sig2 >> sig3 >> sig4;
		
		//Defining cross section
		TF1* XS = new TF1("Cross Section", "[0]*(1+ [1]*x + [2]*(0.5*(3*std::pow(x,2) - 1) ) + [3]*(0.5*(5*std::pow(x,3) - 3*x) ) + [4]*(0.125*(35*std::pow(x,4)-30*x*x+3) ))",-1,1);
		XS->SetParameters(C0, C1, C2, C3, C4, 1, 2, 3, 4);
		
		if (Egamma == E){
			L = XS;
		}
		
	
	}
	
	float Nlam = abs( L->Integral(CthCM(Egamma,pMin), CthCM(Egamma,pMax)) );
	float pathlength = Path_Length(Egamma, L, pMin, pMax);
	
	float Luminosity = Nlam*pathlength*1E-30;  //cm^2;
	
	file.close();
	
	//cout << pathlength << "\t" << Nlam << endl;
	
	if (TMath::IsNaN(Nlam)){
		return 0;
	}
	else{ return Luminosity; }
	
	
	return Luminosity;
}


float Path_Length(float egam, TF1* xs, float lowerMomentum, float upperMomentum){
	int totalLambda = 1000;
	float step = 0.001;
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
		
		//if (TMath::IsNaN(CthL(egam, p))) { cout << "Egamma, plam[low,high]: " << egam << ", [" << lowerMomentum << "," << upperMomentum <<"]" << endl; }
		
		//float cross_section = xs->Eval(CthCM_Kaon(egam, p));
		float theta = acos( CthL(egam, p) );
		float upperLimit = CthCM(egam, upperMomentum);
		float lowerLimit = CthCM(egam, lowerMomentum);

		if (xs_check <= cross_section){
	
			float start = ((float)rand()/(float)(RAND_MAX)) * targetLength;
	
			//float r = ((float)rand()/(float)(RAND_MAX)) *.01;
			float r = ((float)rand()/(float)(RAND_MAX)) *.02 + (-.01);

			//float phi = (((float)rand()/(float)(RAND_MAX)) *2*pi);
			float phi = (((float)rand()/(float)(RAND_MAX)) *2*pi + (-pi));
			
			//float phi_start = (((float)rand()/(float)(RAND_MAX)) *2*pi);
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

			dx = step*cos(phi)*sin(theta);
			dy = step*sin(phi)*sin(theta);
			dz = step*cos(theta);
			while (survivalValue <= prob){
				
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

/********************  TESTED EQUATIONS  ********************/


float CthCM(float egam, float p){
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt(pow(s-xmpr*xmpr,2)/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));

	float Elam = sqrt(p*p+mlam*mlam);
	float Ek = egam + xmpr - Elam;
	float pk = sqrt(Ek*Ek-xmk*xmk);
	
	//float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	float t = -2*xmpr*Elam + xmpr*xmpr + mlam*mlam;
	float cthK = (t+ 2*egam*Ek - xmk*xmk) / (2*egam*abs(pk));
	float angK = acos(cthK);
	//float cthcm = cos( asin((pk/pcm)*sin(angK)) );
	float cthcm = cos( asin((pk/pcm_prime)*sin(angK)) );
	
	if (p > p0(egam)){ cthcm=-cthcm; }
	//if (cthcm<-1.0){ cthcm=-1.0; }
	//if (cthcm>1.0){ cthcm=1.0; }
	
	if (TMath::IsNaN(cthcm)) { cthcm = -1.0; }
	
	return cthcm;
}




float CthL(float egam, float p){	
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	
	float Elam = sqrt(p*p+mlam*mlam);
	float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	float cthcm = CthCM(egam, p);
	//float cthL = -cos( asin((pcm/p)*sin(acos(cthcm))) );
	float cthL = cos( asin((pcm/p)*sin(pi-acos(cthcm))) );
	
	if (p > p0(egam)){ cthL=-cthL; }
	//if (cthL<-1.0){ cthL=-1.0; }
	//if (cthL>1.0){ cthL=1.0; }
	
	return cthL;
}


/************************************************************/

/*
//Calculation partly from PDG
float CthCM(float egam, float p){
	float Elam = sqrt(p*p+mlam*mlam);
	float Ek = egam + xmpr - Elam;
	float pk = sqrt(Ek*Ek-xmk*xmk);

	float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	float s = 2*xmpr*egam+xmpr*xmpr;
	
	float pcm = sqrt(pow(s-xmpr*xmpr,2)/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	
	//Ecm from pdg
	float Ecm = sqrt(xmpr*xmpr + 2*egam*xmpr);
	float xi = log( (pcm+sqrt(xmpr*xmpr+pcm*pcm)/xmpr) );
	
	float t0 = pow( (-xmk*xmk-xmpr*xmpr+mlam*mlam)/(2.0*sqrt(s)) ,2) - pow( pcm-pcm_prime ,2);
	
	float theta = 2.0*asin( sqrt((-t+t0)/ (4.0*pcm*pcm_prime)) );
	
	float cthcm = cos(theta);
	//if (p > p0(egam)){ cthcm=-cthcm; }
	return cthcm;	
	
}
*/

/*
float CthCM(float egam, float p){
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt(pow(s-xmpr*xmpr,2)/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));

	float Elam = sqrt(p*p+mlam*mlam);
	float Ek = egam + xmpr - Elam;
	float pk = sqrt(Ek*Ek-xmk*xmk);
	
	//float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	float t = -2*xmpr*Elam + xmpr*xmpr + mlam*mlam;
	float cthK = (t+ 2*egam*Ek - xmk*xmk) / (2*egam*abs(pk));
	float angK = acos(cthK);
	//float cthcm = cos( asin((pk/pcm)*sin(angK)) );
	float cthcm = cos( asin((pk/pcm_prime)*sin(angK)) );
	
	//if (p > p0(egam)){ cthcm=-cthcm; }
	return cthcm;
}
*/

/*
float CthCM(float egam, float p){	
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt((s-xmpr*xmpr)*(s-xmpr*xmpr)/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	float xi = log( (pcm+sqrt(xmpr*xmpr+pcm*pcm)/xmpr) );
	
	float Elam = sqrt(p*p+mlam*mlam);
	float Ek = egam + xmpr - Elam;
	float pk = sqrt(Ek*Ek-xmk*xmk);
	
	float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	//float cthK = (t+ 2*egam*Ek - xmk*xmk) / (2*egam*abs(pk));
	
	//float cthcm = (pk*cthK - sqrt(pcm*pcm+xmk*xmk)*sinh(xi))/(pcm*cosh(xi));
	
	//float cthcm = ( sqrt(pcm_prime*pcm_prime+mlam*mlam)*cosh(xi)-Elam )/( pcm_prime*sinh(xi) );
	float cthcm = ( sqrt(pcm_prime*pcm_prime+xmk*xmk)*cosh(xi)-Ek )/( pcm_prime*sinh(xi) );
	
	if (p > p0(egam)){ cthcm=-cthcm; }
	//if (cthcm<-1.0){ cthcm=-1.0; }
	//if (cthcm>1.0){ cthcm=1.0; }
	
	return cthcm;
}
*/
/*
float CthCM(float egam, float p){	
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt((s-xmpr*xmpr)*(s-xmpr*xmpr)/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	float xi = log( (pcm+sqrt(xmpr*xmpr+pcm*pcm))/xmpr );
	
	float Elam = sqrt(p*p+mlam*mlam);
	float Ek = egam + xmpr - Elam;
	float pk = sqrt(Ek*Ek-xmk*xmk);
	
	float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	//float cthK = (t+ 2*egam*Ek - xmk*xmk) / (2*egam*abs(pk));
	
	//float cthcm = (pk*cthK - sqrt(pcm*pcm+xmk*xmk)*sinh(xi))/(pcm*cosh(xi));
	
	float cthcm = ( sqrt(pcm_prime*pcm_prime+mlam*mlam)*cosh(xi)-Elam )/( pcm_prime*sinh(xi) );
	//float cthcm = ( sqrt(pcm_prime*pcm_prime+xmk*xmk)*cosh(xi)-Ek )/( pcm_prime*sinh(xi) );
	
	//if (p > p0(egam)){ cthcm=-cthcm; }
	//if (cthcm<-1.0){ cthcm=-1.0; }
	//if (cthcm>1.0){ cthcm=1.0; }
	
	return cthcm;
}
*/
/*
float CthL(float egam, float p){
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt((s-xmpr*xmpr)*(s-xmpr*xmpr)/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	float xi = log( (pcm+sqrt(xmpr*xmpr+pcm*pcm))/xmpr );
	
	float Elam = sqrt(p*p+mlam*mlam);
	
	float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	float cthcm = CthCM(egam, p);
	
	//float cthL = cos( asin((pcm_prime/p)*sin(acos(cthcm))) );
	float cthL = ( -pcm_prime*cthcm*cosh(xi)+sqrt(pcm_prime*pcm_prime+mlam*mlam)*sinh(xi) )/p;
	
	//if (p > p0(egam)){ cthL=-cthL; }
	//if (cthL<-1.0){ cthL=-1.0; }
	//if (cthL>1.0){ cthL=1.0; }
	
	return cthL;
}
*/

/*
float CthL(float egam, float p){
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt((s-xmpr*xmpr)*(s-xmpr*xmpr)/(4*s));
	float pcm_prime = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	float xi = log( (pcm+sqrt(xmpr*xmpr+pcm*pcm))/xmpr );
	
	float Elam = sqrt(p*p+mlam*mlam);
	
	float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	float cthcm = CthCM(egam, p);
	
	float cthL = cos( asin((pcm_prime/p)*sin(acos(cthcm))) );
	
	//if (p > p0(egam)){ cthL=-cthL; }
	//if (cthL<-1.0){ cthL=-1.0; }
	//if (cthL>1.0){ cthL=1.0; }
	
	return cthL;
}
*/

float EgamMin(Float_t P_lambda){
	//Float_t Emin = 0.890938 *(P_lambda) + 0.168662;
	Float_t Emin = (sqrt(P_lambda*P_lambda+mlam*mlam)+xmk-xmpr);
	return Emin;
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
	
	//float Elam = sqrt(pcm_prime*pcm_prime +mlam*mlam)*cosh(xi) - pcm_prime*cthcm*sinh(xi);  //<--When cthcm = 0, what is p?
	float Elam = sqrt(pcm_prime*pcm_prime +mlam*mlam)*cosh(xi);
	float p = sqrt(Elam*Elam-mlam*mlam);

	return p;
}

