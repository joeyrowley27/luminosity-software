
using namespace std;
float pi = TMath::Pi();
float cT = .0789;
float mlam = 1.115;
float xmpr = 0.93827;
float xmk = 0.493677;

float Elam, Ek, pk;
float cthK, cthL, keL, pL, cthcm, angK;

float lum_calc(float Egamma, float pMin, float pMax);

float CthCM(float egam, float p);

float CthL(float egam, float p);

float CthK(float egam, float t);

float EgamMin(Float_t P_lambda);

float Emax(Float_t P_lambda);

float Pmax(Float_t Egam);

float p0(float egam);

void nlamCalc(){

	ofstream outfile;
	outfile.open("Nlam.txt");
	
	int binNumber = 69;
	float Nlum = 0;
	float min, max;
	float Egam_list[]= {1.194, 1.219, 1.244, 1.269, 1.294, 1.32, 1.345, 1.37, 1.395, 1.421, 1.446, 1.471, 1.496, 1.521, 1.546, 1.571, 1.596, 1.621, 1.647, 1.672, 1.698, 1.723, 1.748, 1.774, 1.799, 1.824, 1.849, 1.875, 1.9, 1.925, 1.95, 1.975, 2.001, 2.026, 2.051, 2.076, 2.101, 2.126, 2.151, 2.177, 2.202, 2.227, 2.253, 2.278, 2.303, 2.328, 2.353};
	
	float gflux_list[] = {6.48664046958e+13, 6.48664046958e+13, 6.48664046958e+13, 6.48664046958e+13, 6.48664046958e+13, 6.48664046958e+13, 5.9344693648e+13, 5.35607461185e+13, 4.72070287468e+13, 4.09963504973e+13, 3.55755951537e+13, 3.13535211792e+13, 1.80158320099e+13, 1.52014246255e+13, 1.24939406196e+13, 1.01160708894e+13};
	
	float gflux;
	
	for (int momBin = 0; momBin < 16; momBin++){
		min = 0.6 + 0.1*momBin;
		max = 0.7 + 0.1*momBin;
		
		int startingE;
		for (int item = 1; item <binNumber; item++){
			
			if (Egam_list[item] >= EgamMin(min)){
				startingE = item-1;
				break;
			}
		}
		
		int count = 0;
		float lc;
		float pl =0;
		for(int eBin=startingE; eBin<=46; eBin++){

			lc = lum_calc(Egam_list[eBin], min, max);
			//cout << lc << endl;
			
			if (lc != 0.0){
				Nlum += lc;
				count++;
			}
		}

		gflux = gflux_list[momBin];
		
		outfile << (Nlum/count) << endl;
		
		Nlum = 0;
		
	}

	outfile.close();
}


float lum_calc(float Egamma, float pMin, float pMax){
	
	ifstream file ("KLambda_legendre.txt");
	float s = 2*xmpr*Egamma+xmpr*xmpr;
	float pcm = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	
	gSystem->Load("libMathMore");
	string line;
	float Index, E, W, Order, C0, C1, C2, C3, C4, sig0, sig1, sig2, sig3, sig4;
	
	
	TF1* L;
	int count_line = 0;
	while (!file.eof()){
		getline(file,line);
		count_line++;
		if (count_line >= 69 && count_line <= 150){

			file >> Index >> E >> W >> Order >> C0 >> C1 >> C2 >> C3 >> C4 >> sig0 >> sig1 >> sig2 >> sig3 >> sig4;
			
			//Defining cross section
			TF1* XS = new TF1("Cross Section", "[0]*(1+[1]*ROOT::Math::legendre([5],x)+[2]*ROOT::Math::legendre([6],x)+[3]*ROOT::Math::legendre([7],x)+[4]*ROOT::Math::legendre([8],x))", -1, 1);
			XS->SetParameters(C0, C1, C2, C3, C4, 1, 2, 3, 4);
			
			if (Egamma == E){
				L = XS;
			}
		}
	
	}
	
	float Nlam = abs( L->Integral(CthCM(Egamma,pMin), CthCM(Egamma,pMax)) )*(6E23*0.071*40.0)*1E-30;  //cm^2
	
	if (TMath::IsNaN(Nlam)){
		return 0;
		
	}
	
	file.close();
	
	return Nlam;
}



float CthL(float egam, float p){
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt(((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s));
	
	float t = -(2*xmpr*Elam - xmpr*xmpr - mlam*mlam);
	Elam = sqrt(p*p+mlam*mlam);
	pL = sqrt(Elam*Elam-mlam*mlam);
	float cthcm = CthCM(egam, p);
	cthL = cos( asin((pcm/pL)*sin(acos(cthcm))) );
	
	return cthL;
}

float CthK(float egam, float t){
	t = -abs(t);
	Elam = (xmpr*xmpr + mlam*mlam - t) / (2*xmpr);
	Ek = egam + xmpr - Elam;
	pk = sqrt(Ek*Ek-xmk*xmk);
	cthK = (t+ 2*egam*Ek - xmk*xmk) / (2*egam*abs(pk));
	return cthK;
}

float CthCM(float egam, float p){
	float s = 2*xmpr*egam+xmpr*xmpr;
	float pcm = sqrt((s-xmpr*xmpr)*(s-xmpr*xmpr)/(4*s));
	float pcm_prime = sqrt( ((s-xmk*xmk-mlam*mlam)*(s-xmk*xmk-mlam*mlam)-4*xmk*xmk*mlam*mlam)/(4*s) );
	float xi = log( (pcm+sqrt(xmpr*xmpr+pcm*pcm))/xmpr );
	
	float Elam = sqrt(p*p+mlam*mlam);
	float Ek = egam + xmpr - Elam;
	float pk = sqrt(Ek*Ek-xmk*xmk);
	
	float cthcm = ( sqrt(pcm_prime*pcm_prime+xmk*xmk)*cosh(xi)-Ek )/( pcm_prime*sinh(xi) );
	
	if (p > p0(egam)){ cthcm=-cthcm; }
	if (cthcm<-1.0){ cthcm=-1.0; }
	if (cthcm>1.0){ cthcm=1.0; }
	
	return cthcm;
}


float EgamMin(Float_t P_lambda){
	Float_t Emin = 0.890938 *(P_lambda) + 0.168662;
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

