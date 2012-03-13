
void FillAsnAntiPar(char *Asn1, char *Asn2, CHAIN **Chain, int Cn1, int Cn2, 
		    PATTERN **Pat, int NPat, StrideCmd *Cmd);
void FillAsnPar(char *Asn1, char *Asn2, CHAIN **Chain, int Cn1, int Cn2, 
		PATTERN **Pat, int NPat, StrideCmd *Cmd);
void FilterAntiPar(PATTERN **Pat, int NPat);
void FilterPar(PATTERN **Pat, int NPat);
