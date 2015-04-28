#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gl/glut.h>
#include <gl/gl.h>
#include <gl/glu.h>
#include <gl/glui.h>
#include <vector>
#include <iostream>
#include <tchar.h>
#include <sstream>

//nasconde la console
#pragma comment(linker,"/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

//definizione funzioni base di Hermite
#define PHI0(t) (2.0*t*t*t - 3.0*t*t + 1)
#define PHI1(t) (t*t*t - 2.0*t*t + t)
#define PSI0(t) (-2.0*t*t*t + 3.0*t*t)
#define PSI1(t) (t*t*t - t*t)

#define PESOBASE 1.0

int larghezza_principale = 900, altezza_principale = 800;
int larghezza_secondaria = 400, altezza_secondaria = 300;
int punto_selezionato = -1;

int winIdPrincipale, winIdFunzioniBase; //ID Finestre

int metodo_attivo = 0; //1: Hermite; 2: Bezier; 3: curve Spline
int modifica_derivata = 0; //0: Non modifica derivate, 1: Modifica derivate
int scelta_opzioni = 0, modifica_molteplicità = 0, indice_nodo, valore_molteplicità;
int ordine_Spline = 4; //Ordine spline

float valore_suddivisione; //valore di suddivisione
int attiva_suddivisione = 0; //1: viene applicato l'algoritmo subdivision a bezier

GLUI_Panel *pannello_opzioni;
GLUI_RadioGroup *radio_opzioni;

//pannello scelta parametri
int scelta_parametrizzazione = 0;
GLUI_Panel *pannello_parametrizzazione;
GLUI_RadioGroup *radio_parametrizzazione;

//pannello per esecuzione Hermite
GLUI_Panel *pannello_Hermite;
GLUI_Button *bottone_Hermite;

//pannello per esecuzione Bezier
GLUI_Panel *pannello_Bezier;
GLUI_Button *buttonBezier;

//pannello per esecuzione Spline
int scelta_parametrizzazione_spline = 0;
GLUI_Panel *pannello_Spline;
GLUI_Button *bottone_Spline;

//spinner subdivision
GLUI_Spinner *spinner_subd;

//spinner molteplicità
GLUI_Spinner *spinner_i_nodo, *spinner_molteplicità;

using namespace std;

typedef struct glPoint2D{
	GLfloat x,y;
}GLPOINT2D;

//contenitore di dati
vector <GLPOINT2D> Punti;

//contenitore dei pesi dei punti
vector <float> PesiPunti;

//conserva le coordinate dei punti dove è stata modificata la derivata
vector <GLPOINT2D> DerivateMod;

//Punti preparati per controllare se va tutto bene
vector <GLPOINT2D> PuntiPrepatati;

//converte un intero in stringa
string int2str(int x) 
{
	stringstream ss;
	ss << x;
	return ss.str();
}

//per acquisire l'input del mouse
void myMouse(int button, int state, GLint xmouse, GLint ymouse){

	int tolleranzaDistanzaClick = 3;
	float distanza, distanza1;

	GLPOINT2D newPoint, zero;
	zero.x = 0.0;
	zero.y = 0.0;
	newPoint.x = xmouse;
	newPoint.y = altezza_principale - ymouse;

	if (state == GLUT_DOWN){ //se lo stato del bottone è premuto
		switch(button){

		case GLUT_LEFT_BUTTON:

			//scelta_opzioni = 0 --> modalita' inserimento
			//scelta_opzioni = 1 --> modifica punti inseriti
			//scelta_opzioni = 2 --> modifica pesi
			//modifica_derivata = 0 --> se non voglio modificare la derivata sui punti inseriti
			//modifica_derivata = 1 --> se voglio modificare la derivata sui punti inseriti
			if(scelta_opzioni == 1 || modifica_derivata == 1 || scelta_opzioni == 2){

				//ho già inserito dei punti, devo "catturare" il punto cliccato
				//tramite un controllo "ad area" che confronta la zona cliccata
				//coi punti inseriti, se la distanza è minore di una costante
				//allora selezionerò il punto per poi modificarne la posizione
				if(Punti.size() > 0){
					punto_selezionato = 0;

					//calcoliamo la distanza dal newpoint da tutti i punti inseriti prima e prendo l'indice del punto più vicino
					distanza = sqrt((Punti.at(0).x - newPoint.x)*(Punti.at(0).x - newPoint.x) + (Punti.at(0).y - newPoint.y)*(Punti.at(0).y - newPoint.y));

					for (int i = 1; i < Punti.size(); i++){
						distanza1 = sqrt((Punti.at(i).x - newPoint.x)*(Punti.at(i).x - newPoint.x) + (Punti.at(i).y - newPoint.y)*(Punti.at(i).y - newPoint.y));

						//faccio il controllo sui minimi delle distanze cioè prendo quello con distanza minore
						if (distanza1 < distanza){
							punto_selezionato = i;
							distanza = distanza1;
						}
					}

					//confronto ora col mio indice di tolleranza
					if(distanza > tolleranzaDistanzaClick){
						punto_selezionato = -1;
					}
				}
			} else if (scelta_opzioni == 0){

				//inserisce il nuovo punto ed i dati relativi a quel punto
				Punti.push_back(newPoint);
				PesiPunti.push_back(PESOBASE);
				DerivateMod.push_back(zero);
			}
			break;

			//tasto destro mouse elimino tutti i punti inseriti
		case GLUT_RIGHT_BUTTON:
			Punti.clear();
			PesiPunti.clear();
			DerivateMod.clear();
			break;

			//Tasto centrale del mouse elimino l'ultimo punti inserito
		case GLUT_MIDDLE_BUTTON:
			if (Punti.size() > 0){
				Punti.pop_back();
				PesiPunti.pop_back();
				DerivateMod.pop_back();
			}
			break;
		}
	} else {
		//lo stato del bottone non è premuto 
		switch(button){
		case GLUT_LEFT_BUTTON:
			punto_selezionato = -1;
			break;
		}

	}
	glutPostRedisplay();
}

//utilizzato per spostare un punto oppure per modificare la sua derivata
void mouseMove(GLint xmouse, GLint ymouse){

	GLPOINT2D newPoint;
	newPoint.x = xmouse;
	newPoint.y = altezza_principale - ymouse;
	
	if(punto_selezionato >= 0){
		if(scelta_opzioni == 1){
			//sostituisco il punto selezionato con quello nuovo
			Punti.at(punto_selezionato) = newPoint;
		}else if(scelta_opzioni == 2 && (metodo_attivo == 2 || metodo_attivo == 3)){
			//controllo la coordinata y del punto selezionato e del punto corrente per modificare il peso del punto selezionato
			if (newPoint.y > Punti.at(punto_selezionato).y){
				PesiPunti.at(punto_selezionato) = (newPoint.y - Punti.at(punto_selezionato).y) / 15;
			}else if (newPoint.y < Punti.at(punto_selezionato).y){
				if (PesiPunti.at(punto_selezionato) - (Punti.at(punto_selezionato).y - newPoint.y) <= 0){
					PesiPunti.at(punto_selezionato) = PESOBASE - 1;
				}else{
					PesiPunti.at(punto_selezionato) = (Punti.at(punto_selezionato).y - newPoint.y) / 15;
				}
			}
		} else if (modifica_derivata == 1){

			//la moltiplicazione per 5 è una regola
			float derx = (newPoint.x - Punti.at(punto_selezionato).x)*5;
			float dery = (newPoint.y - Punti.at(punto_selezionato).y)*5;

			DerivateMod.at(punto_selezionato).x = derx;
			DerivateMod.at(punto_selezionato).y = dery;
		}
	}
	glutPostRedisplay();
}

//genera t (che punta ad una variabile creata in display)
void parametrizzazione_uniforme(float* t){

	//definisco il passo della parametrizzazione, dividendo l'ampiezza dell'intervallo
	//per il numero di punti presenti -1 (l'ampiezza dell'intervallo è sempre 1)
	float step = 1.0/(float)(Punti.size()-1);

	//moltiplica l'indice per il passo
	//es. con 5 punti => t[0] = 0, t[1] = 0.25, t[2] = 0.5, t[3] = 0.75, t[4] = 1
	for (int i = 0; i < Punti.size(); i++){
		t[i] = i * step;
	}

}

//genera t (che punta ad una variabile creata in display)
void parametrizzazione_corde(float* t){

	t[0] = 0;

	//prende il valore t precedente e lo somma alla distanza tra il punto corrente del ciclo e quello precedente (lunghezza del segmento che li unisce)
	for (int i = 1; i<Punti.size(); i++){
		t[i] = t[i-1] + sqrt((Punti.at(i).x - Punti.at(i-1).x)*(Punti.at(i).x - Punti.at(i-1).x) + (Punti.at(i).y - Punti.at(i-1).y)*(Punti.at(i).y - Punti.at(i-1).y));
	}

	//divido per il massimo dell'ultimo componente per riportarlo nel range 0-1 come nella parametrizzazone_uniforme (ma valori diversi)
	for (int i = 0; i<Punti.size(); i++){
		t[i] = t[i]/t[Punti.size()-1];
	}

}

//derivate rispetto a t della componente parametrica in x
//metodo del rapporto incrementale
float dx(int i, float* t){
	//i??
	if (i <= 0){
		return 0;
	}

	//altrimenti restituisco la derivata in x con il metodo del rapporto incrementale
	return (Punti.at(i).x - Punti.at(i-1).x)/(t[i] - t[i-1]);
}

//derivate rispetto a t della componente parametrica in y
//metodo del rapporto incrementale
float dy(int i, float* t){
	if (i <= 0){
		return 0;
	}

	//altrimenti restituisco la derivata in y con il metodo del rapporto incrementale
	return (Punti.at(i).y - Punti.at(i-1).y)/(t[i] - t[i-1]);
}

//utilizzato dentro interpolazione_Hermite
float DX(int i, float *t){

	//non è stata calcolata la derivata del punto i
	if(DerivateMod.at(i).x == 0){
		return dx(i, t);
	}

	//è già stata calcolata la derivata del punto i
	if(DerivateMod.at(i).x != 0){
		return DerivateMod.at(i).x;
	}
}

//utilizzato dentro interpolazione_Hermite
float DY(int i, float *t){

	//non è stata calcolata la derivata del punto i
	if(DerivateMod.at(i).y == 0){
		return dy(i, t);
	}

	//è già stata calcolata la derivata del punto i
	if(DerivateMod.at(i).y != 0){
		return DerivateMod.at(i).y;
	}
}

//implementazione interpolazione di Hermite
void interpolazione_Hermite(float* t){

	//valutiamo la nostra curva su 1000 valori
	//numero di valori del parametro t in cui valutare la curva interpolante di Hermite
	int nvpt = 1000;

	float passot = 1.0/(float)(nvpt-1); //0,001001001001001

	//per ogni valore, devo vedere in quale sottointervallo cade il mio punto di valutazione
	int is = 0; //Indice del sottointervallo a cui appartiene il valore del parametro in cui valutare la curva
	glColor3f(0.0,0.4,0.0);
	glBegin(GL_LINE_STRIP);
	//1000 valori tra i due estremi che sarebbero i due punti tra i quali viene disegnata
	for (float ti = 0; ti <= 1; ti += passot){

		if (ti > t[is+1]){
			is++; //is è l'indice dell'intervallo a cui appartiene il valore ti del parametro
		}

		//mappare il valore ti (appartenente all'intervallo [t[is],t[is+1]) nell'intervallo [0,1]
		//perchè le funzioni base sono state definite su un intervallo [0,1]
		float tim; //tim è la trasformazione affine che apparterrà sempre all'intervallo [0,1]
		float range = (t[is+1] - t[is]);
		tim = (ti - t[is])/range; 

		GLPOINT2D IH;
		
		//PHI0, PHI1, PSI0, PSI1 sono le funzioni base di Hermite
		IH.x = Punti.at(is).x*PHI0(tim) + DX(is,t)*PHI1(tim)*range + Punti.at(is+1).x*PSI0(tim) + DX(is+1,t)*PSI1(tim)*range;
		IH.y = Punti.at(is).y*PHI0(tim) + DY(is,t)*PHI1(tim)*range + Punti.at(is+1).y*PSI0(tim) + DY(is+1,t)*PSI1(tim)*range;

		//disegna uno dei mille valori interni ai due estremi
		glVertex2f(IH.x, IH.y);
	}
	glEnd();
}


void funzione_base_Hermite(){

	glColor3f(0.0, 0.0, 0.0);
	glRasterPos2f(0.42, 0.97);
	string titolo = "Hermite";
	int len = titolo.length();
	for (int l = 0; l < len; l++){
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, titolo[l]);
	}

	//valutiamo la nostra curva su 1000 valori
	//numero di valori del parametro t in cui valutare la curva interpolante di Hermite
	int nvpt = 1000;

	float passot = 1.0/(float)(nvpt-1);

	//4 volte perchè le funzioni base sono 4
	for (int i = 0; i < 4;i++){
		glBegin(GL_LINE_STRIP);

		//la coordinata x sarà sempre il valore ti invece la coordinata y sarà calcolata in base alla funzione base (sono 4) presa in considerazione
		for (float ti = 0; ti <= 1; ti += passot){
			GLPOINT2D IH;
			IH.x = ti;

			if (i==0){

				glColor3f(1.0,0.0,0.0);
				IH.y = PHI0(ti);

			}else if(i == 1){

				glColor3f(0.0,0.4,0.0);
				IH.y = PHI1(ti);

			}else if(i == 2){

				glColor3f(0.0,0.0,1.0);
				IH.y = PSI0(ti);

			}else if(i == 3){

				glColor3f(0.0,0.0,0.0);
				IH.y = PSI1(ti);
			}
			glVertex2f(IH.x, IH.y);
		}
		glEnd();
	}
}

//algoritmo di suddivisione
void Subdivision(){

	GLPOINT2D *c = new GLPOINT2D[Punti.size()];
	GLPOINT2D *c1 = new GLPOINT2D[Punti.size()];
	GLPOINT2D *c2 = new GLPOINT2D[Punti.size()];

	for(int i = 0; i < Punti.size() ; i++)
		c[i] = Punti.at(i);

	c1[0] = c[0];
	c2[0] = c[Punti.size()-1];

	for(int j = 1; j <= Punti.size(); j++){
		for(int i = 0; i < Punti.size() - j ; i++){
			c[i].x = c[i].x * (1 - valore_suddivisione) + valore_suddivisione * c[i+1].x;
			c[i].y = c[i].y * (1 - valore_suddivisione) + valore_suddivisione * c[i+1].y;
		}

		c1[j].x = c[0].x;
		c1[j].y = c[0].y;
		c2[j] = c[Punti.size() - 1 - j];
	}

	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for(int i = 0; i < Punti.size() ; i++)
		glVertex2f(c1[i].x, c1[i].y);
	glEnd();

	glBegin(GL_LINE_STRIP);
	for(int i = 0; i < Punti.size() ; i++)
		glVertex2f(c2[i].x, c2[i].y);
	glEnd();

}

//implementazione interpolazione di Bezier
void Bezier(){

	//valutiamo la nostra curva su 1000 valori
	int npv = 1000;
	float tstep = 1.0/(float)(npv - 1);

	//creati due vettori della stessa dimensione dei punti
	GLPOINT2D *c = new GLPOINT2D[Punti.size()];
	float *w = new float[PesiPunti.size()];

	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 0.0, 0.0);
	for(float ti = 0; ti <= 1; ti += tstep){

		//inserisce nel vettore c la moltiplicazione tra il peso per le coordinate x e y di tutti i punti
		//e nel vettore w inserisce solo il peso dei punti
		for(int i = 0; i < Punti.size() ; i++)
		{
			c[i].x = PesiPunti.at(i) * Punti.at(i).x;
			c[i].y = PesiPunti.at(i) * Punti.at(i).y;
			w[i] = PesiPunti.at(i);
		}

		//utilizza algoritmo De Casteljau applicato ai punti ed ai pesi
		for(int j = 1; j <= Punti.size(); j++){
			for(int i = 0; i < Punti.size() - j ; i++){
				c[i].x = c[i].x * (1 - ti) + ti * c[i+1].x;
				c[i].y = c[i].y * (1 - ti) + ti * c[i+1].y;
				w[i] = w[i] * (1 - ti) + ti * w[i+1];
			}
		}
		//disegna solo il primo dei valori che sarà il punto desiderato
		//per avere la curva in 2D è necessario dividere per la terza componente (peso) perchè quando il peso è maggiore di 1 abbiamo una curva razionale (3D)
		glVertex2f(c[0].x / w[0], c[0].y / w[0]);
	}
	glEnd();

	delete(c);
	delete(w);
}

//resta da capire cosa sia distanza 1 beta (d1b)
void funzione_base_Bezier(){

	glColor3f(0.0, 0.0, 0.0);
	glRasterPos2f(0.44, 0.97);
	string titolo = "Bezier";
	int len = titolo.length();
	for (int l = 0; l < len; l++){
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, titolo[l]);
	}

	int scelta_colore = 0;
	int campioni = 300;
	float ti = 0;
	float passot = 1.0/(campioni - 1);
	float *sommaPesi = new float[campioni]; 

	//allochiamo la matrice B, nella quale salveremo le nostre funzioni base
	float **B = new float *[campioni];

	//posizione dove memorizzare la valutazione che corrisponde al valore del parametro t
	for (int k = 0; k < campioni; k++, ti += passot){
		B[k] = new float[Punti.size()+1]();
		B[k][Punti.size()] = 1.0; //condizione iniziale

		//formula ottimizzata, per evitare controlli su valori nulli
		//ricorsiva, ci salviamo solo il valore che ci interessa
		for (int i = 1; i < Punti.size(); i++){
			float d1b = 0; //distanza 1 beta
			for (int j = 0; j < Punti.size(); j++){
				B[k][j] = (1-ti) * B[k][j+1] + d1b;
				d1b = ti * B[k][j+1];
			}
			B[k][Punti.size()] = d1b;
		}

		sommaPesi[k] = 0;
		for (int j = 0; j < Punti.size(); j++){
			sommaPesi[k] += B[k][j+1] * PesiPunti.at(j);
		}

	}

	//disegnamo le funzioni base, che sono le colonne della matrice B
	glColor3f(1.0, 0.0, 0.0);
	int cont = 0;
	for (int ibase = 0; ibase < Punti.size(); ibase++){

		switch (scelta_colore)
		{
		case 0:
			scelta_colore++;
			glColor3f(0.0, 1.0, 0.0);
			break;
		case 1:
			scelta_colore++;
			glColor3f(0.0, 0.5, 1.0);
			break;
		case 2:
			scelta_colore++;
			glColor3f(1.0, 0.0, 0.5);
			break;
		case 3:
			scelta_colore++;
			glColor3f(1.0, 0.5, 0.0);
			break;
		case 4:
			scelta_colore++;
			glColor3f(0.6, 0.0, 0.6);
			break;
		case 5:
			scelta_colore = 0;;
			glColor3f(0.2, 1.0, 0.6);
			break;
		}

		glBegin(GL_LINE_STRIP);
		for (int l = 0; l < campioni; l++){
			glVertex2f(passot*(float)l, B[l][ibase+1]*PesiPunti.at(ibase)/sommaPesi[l]);
		}
		glEnd();
	}

	delete(sommaPesi);
	delete(B);
}

//Utilizzata per creare la partizione nodale estesa
void costruisci_nodi(float *t, float *Nodi, char* molteplicità)
{
	int i, cont;
	int k = Punti.size() - ordine_Spline; //Numero di Nodi interni all'intervallo

	//Nodi fittizi a sinistra
	for (i = 0; i < ordine_Spline; i++)
	{
		Nodi[i] = 0;
		molteplicità[i] = '4';
	}

	//Costruzione nodi veri
	cont = 2;
	for (i = ordine_Spline; i < ordine_Spline + k; i++)
	{
		Nodi[i] = t[cont];
		molteplicità[i] = '1';
		cont++;
	}

	//Nodi fittizi a destra
	for (i = ordine_Spline + k; i < 2 * ordine_Spline + k; i++)
	{
		Nodi[i] = 1;
		molteplicità[i] = '4';
	}

	if (modifica_molteplicità == 1 && Punti.size() > 4)
	{
		float val_nodo = Nodi[indice_nodo];
		if (valore_molteplicità == 2)
		{
			molteplicità[indice_nodo] = '2';
		}
		if (valore_molteplicità == 3)
		{
			molteplicità[indice_nodo] = '3';
		}
		if (valore_molteplicità == 4)
		{
			molteplicità[indice_nodo] = '4';
		}
		if (valore_molteplicità > 1)
		{
			for (i = 1; i < valore_molteplicità; i++)
			{
				if (indice_nodo + i == ordine_Spline + k){
					break;
				}else{
					Nodi[indice_nodo + i] = val_nodo;
					molteplicità[indice_nodo + i] = molteplicità[indice_nodo];
				}
			}
		}
	}
}


int localizza_intervallo_internodale(float t, float *Nodi)
{
	//Implementazione del metodo di bisezione
	int a = ordine_Spline - 1;
	int b = Punti.size();

	int index;

	while(b - a > 1)
	{
		index = (a + b) / 2;
		if (t < Nodi[index])
		{
			b = index;
		} else {
			a = index;
		}
	}
	return a;
}

void disegna_base_Spline(float *Nodi, char *molteplicità)
{
	glColor3f(0.0, 0.0, 0.0);
	glRasterPos2f(0.43, 0.97);
	string titolo = "Spline";
	int len = titolo.length();
	for (int l = 0; l < len; l++){
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, titolo[l]);
	}

	int scelta_colore = 0;
	int Ncampioni = 300;
	float ti = 0;
	float dt = 1.0/(float)(Ncampioni - 1);

	float *sommaPesi = new float[Ncampioni];

	float **B = new float*[Punti.size()];

	for (int i = 0; i < Punti.size(); i++) {
		B[i] = new float[Ncampioni]();
	}
	//valutazione di ciascuna funzione base per ogni valore del parametro t mediante le formule di Cocks
	for (int k = 0; k < Ncampioni; k++, ti += dt)
	{
		int l = localizza_intervallo_internodale(ti, Nodi);

		B[l][k] = 1;
		for (int i = 0; i < ordine_Spline - 1; i++)
		{
			float tmp = 0.0;
			for (int j = l - i; j <= l; j++)
			{
				float d1 = ti - Nodi[j];
				float d2 = Nodi[i+j+1] - ti;
				float beta = B[j][k] / (d1 + d2);
				B[j-1][k] = d2 * beta + tmp;
				tmp = d1 * beta;
			}
			B[l][k] = tmp;
		}

		sommaPesi[k] = 0;
		for (int j = 0; j < Punti.size(); j++){
			sommaPesi[k] += B[j][k] * PesiPunti.at(j);
		}

	}

	for (int i = 0; i < Punti.size(); i++)
	{
		
		switch (scelta_colore)
		{
		case 0:
			scelta_colore++;
			glColor3f(0.0, 1.0, 0.0);
			break;
		case 1:
			scelta_colore++;
			glColor3f(0.0, 0.5, 1.0);
			break;
		case 2:
			scelta_colore++;
			glColor3f(1.0, 0.0, 0.5);
			break;
		case 3:
			scelta_colore++;
			glColor3f(1.0, 0.5, 0.0);
			break;
		case 4:
			scelta_colore++;
			glColor3f(0.6, 0.0, 0.6);
			break;
		case 5:
			scelta_colore = 0;;
			glColor3f(0.2, 1.0, 0.6);
			break;
		}

		ti = 0; 
		glBegin(GL_LINE_STRIP);
		for (int k = 0; k < Ncampioni; k++, ti += dt)
		{
			glVertex2f(ti, (B[i][k]*PesiPunti.at(i))/sommaPesi[k]);
		}
		glEnd();
	}

	glColor3f(0.0, 0.4, 0.0);
	for (int j = 0; j < Punti.size() + ordine_Spline; j++)
	{
		glRasterPos2f(Nodi[j]-0.01, -0.05);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, molteplicità[j]);
	}

	delete(sommaPesi);
	delete(B);
}

//Utilizzata per disegnare la spline
void DeBoor(float *t, float *Nodi)
{
	int nvalorit = 1000;
	float tstep = 1.0 / (float)(nvalorit - 1);

	GLPOINT2D *c = new GLPOINT2D[ordine_Spline];
	float *w = new float[ordine_Spline];

	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINE_STRIP);
	for (float vt = 0; vt <= 1; vt += tstep)
	{
		int l = localizza_intervallo_internodale(vt, Nodi);

		//Implementamo l'algoritmo di DeBoor
		for (int i = 0; i < ordine_Spline; i++)
		{
			c[i].x = Punti.at(i + l - ordine_Spline + 1).x * PesiPunti.at(i + l - ordine_Spline + 1);
			c[i].y = Punti.at(i + l - ordine_Spline + 1).y * PesiPunti.at(i + l - ordine_Spline + 1);
			w[i] = PesiPunti.at(i + l - ordine_Spline + 1);
		}
		for (int j = 0; j < ordine_Spline - 1; j++)
		{
			for (int i = ordine_Spline - 1; i > j; i--)
			{
				int ti = l - ordine_Spline + 1 + i;
				int timj = l + i - j;

				float den = Nodi[timj] - Nodi[ti];
				float dt = (vt - Nodi[ti]) / den;

				c[i].x = c[i].x * dt + (c[i-1].x * (1 - dt));
				c[i].y = c[i].y * dt + (c[i-1].y * (1 - dt));
				w[i] = w[i] * dt + (w[i-1] * (1 - dt));
			}
		}
		glVertex2f(c[ordine_Spline - 1].x/w[ordine_Spline - 1],c[ordine_Spline - 1].y / w[ordine_Spline - 1]);
	}
	glEnd();

	delete(w);
	delete(c);
}

//posiziona punti preparati
void inizializza_punti(){

	Punti.clear();
	PesiPunti.clear();
	DerivateMod.clear();

	GLPOINT2D zero = {0.0,0.0};
	for (int i = 0; i < PuntiPrepatati.size(); i++){
		Punti.push_back(PuntiPrepatati.at(i));
		PesiPunti.push_back(PESOBASE);
		DerivateMod.push_back(zero);
	}
}

void scelta_metodi(int scelta){

	switch(scelta){
	case 1:
		metodo_attivo = 1; //Hermite
		glutPostRedisplay();
		break;

	case 2:
		metodo_attivo = 2;	//Bezier
		glutPostRedisplay();
		break;

	case 3:
		metodo_attivo = 3; //Spline
		glutPostRedisplay();
		break;

	case 4:
		metodo_attivo = 0; //nessun metodo attivo
		inizializza_punti(); //inserisco punti già preparati
		glutPostRedisplay();
		break;
	}

}

void display(){

	glutSetWindow(winIdFunzioniBase);
	glClear(GL_COLOR_BUFFER_BIT); //pulisco la finestra secondaria (delle funzioni base)

	glutSetWindow(winIdPrincipale);
	glClear(GL_COLOR_BUFFER_BIT); //pulisco la finestra principale (dei punti)

	//modifica il limite dello spinner per la scelta dei nodi
	if (Punti.size() > ordine_Spline){
		spinner_i_nodo -> set_int_limits(ordine_Spline,Punti.size()-1);
	}else{
		spinner_i_nodo -> set_int_limits(0,0);
	}

	//disegno i punti inseriti fin'ora
	glBegin(GL_POINTS);
	for (int i = 0; i < Punti.size(); i++){
		if (punto_selezionato == i){
			glColor3f(0.0, 0.0, 1.0); 
		}else{
			glColor3f(0.0, 0.0, 0.0); 
		}
		glVertex2f(Punti.at(i).x, Punti.at(i).y);
	}
	glEnd();

	//stampo gli indici dei punti
	for (int i = 0; i < Punti.size(); i++){
		if (punto_selezionato == i){
			glColor3f(0.0, 0.0, 1.0); 
		}else{
			glColor3f(0.0, 0.0, 0.0); 
		}
		glRasterPos2f(Punti.at(i).x - 15, Punti.at(i).y + 10);
		string indice = int2str(i+1);
		int len = indice.length();
		for (int l = 0; l < len; l++){
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, indice[l]);
		}
	}

	glColor3f(1.0,0.0,0.0);
	glLineStipple(5, 0xAAAA);
	glEnable(GL_LINE_STIPPLE);
	glBegin(GL_LINE_STRIP);	//disegno le linee che collegano i punti
	for (int i = 0; i < Punti.size(); i++)
		glVertex2f(Punti.at(i).x, Punti.at(i).y);
	glEnd();
	glDisable(GL_LINE_STIPPLE);

	//Parametrizzazione per 
	float* t = new float[Punti.size()];

	if (scelta_parametrizzazione == 0){
		parametrizzazione_uniforme(t);
	} else if (scelta_parametrizzazione == 1){
		parametrizzazione_corde(t);
	}

	if(Punti.size() > 1){

		//curve interpolanti di Hermite
		if (metodo_attivo == 1){

			interpolazione_Hermite(t);

			if(punto_selezionato >= 0 && modifica_derivata == 1){
				//il punto selezionato
				GLPOINT2D P0 = Punti.at(punto_selezionato);

				//coordinate del punto corrispondente al valore del parametro t = 0
				//appartenente alla retta tangente alla curva nel punto selezionato
				GLPOINT2D P1;

				P1.x = P0.x - DerivateMod.at(punto_selezionato).x * t[punto_selezionato];
				P1.y = P0.y - DerivateMod.at(punto_selezionato).y * t[punto_selezionato];

				//coordinate del punto corrispondente al valore del parametro t = 1
				//appartenente alla retta tangente alla curva nel punto selezionato
				GLPOINT2D P2;

				P2.x = P0.x + DerivateMod.at(punto_selezionato).x - DerivateMod.at(punto_selezionato).x * t[punto_selezionato];
				P2.y = P0.y + DerivateMod.at(punto_selezionato).y - DerivateMod.at(punto_selezionato).y * t[punto_selezionato];

				glBegin(GL_LINE_STRIP);
				glVertex2f(P1.x, P1.y);
				glVertex2f(P2.x, P2.y);
				glEnd();
			}

			glutSetWindow(winIdFunzioniBase);
			funzione_base_Hermite();

		} else if (metodo_attivo == 2){

			Bezier();
			if(attiva_suddivisione == 1){
				Subdivision();
			}
			glutSetWindow(winIdFunzioniBase);
			funzione_base_Bezier();

		} else if (metodo_attivo == 3) {
			if (Punti.size() >= ordine_Spline) {

				float *Nodi = new float[Punti.size() + 2 * ordine_Spline];
				char *molt = new char[Punti.size() + ordine_Spline];
				float* t_spline = new float[Punti.size()];
				if (scelta_parametrizzazione_spline == 0){
					parametrizzazione_uniforme(t_spline);
				}else{
					parametrizzazione_corde(t_spline);
				}
				costruisci_nodi(t_spline, Nodi, molt);
				DeBoor(t_spline, Nodi);

				glutSetWindow(winIdFunzioniBase);
				disegna_base_Spline(Nodi, molt);

				delete(t_spline);
				delete(Nodi);
				delete(molt);
			}
		}
	}

	glutSetWindow(winIdPrincipale);
	glutSwapBuffers();
	glutSetWindow(winIdFunzioniBase);
	glutSwapBuffers();
	delete(t);
}

void myinit (void)
{
	glutSetWindow(winIdFunzioniBase);
	glClearColor(1.0, 1.0, 1.0, 0.0); //colore dello sfondo finestra secondaria

	//Definisco il sistema di riferimento per la finestra secondaria di disegno delle funzioni base
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-0.05, 1.05, -0.18, 1.05);
	glutSetWindow(winIdPrincipale);
	glClearColor(1.0, 1.0, 1.0, 0.0); //colore dello sfondo finestra principale
	glPointSize(5.0);
	//Definisco il sistema di riferimento per la finestra principale di disegno delle curve
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,float(larghezza_principale),0.0,float(altezza_principale));

	GLPOINT2D newPoint = {118,352};
	PuntiPrepatati.push_back(newPoint);
	newPoint.x = 241;
	newPoint.y = 544;
	PuntiPrepatati.push_back(newPoint);
	newPoint.x = 487;
	newPoint.y = 543;
	PuntiPrepatati.push_back(newPoint);
	newPoint.x = 457;
	newPoint.y = 342;
	PuntiPrepatati.push_back(newPoint);
	newPoint.x = 258;
	newPoint.y = 231;
	PuntiPrepatati.push_back(newPoint);
	newPoint.x = 422;
	newPoint.y = 143;
	PuntiPrepatati.push_back(newPoint);
	newPoint.x = 657;
	newPoint.y = 231;
	PuntiPrepatati.push_back(newPoint);
	newPoint.x = 521;
	newPoint.y = 285;
	PuntiPrepatati.push_back(newPoint);
}

void createOptionGlui(){

	GLUI *glui = GLUI_Master.create_glui_subwindow(winIdPrincipale, GLUI_SUBWINDOW_TOP);

	//Hermite
	pannello_Hermite = glui->add_panel("Hermite", GLUI_PANEL_EMBOSSED);
	bottone_Hermite = glui->add_button_to_panel(pannello_Hermite,"HERMITE",1,scelta_metodi);
	radio_parametrizzazione = glui->add_radiogroup_to_panel(pannello_Hermite, &scelta_parametrizzazione);
	glui->add_radiobutton_to_group(radio_parametrizzazione, "Parametrizzazione Uniforme");
	glui->add_radiobutton_to_group(radio_parametrizzazione, "Parametrizzazione delle Corde");

	//possibilità di modificare le derivate per modificare il comportamento delle nostre curve interpolate
	glui->add_checkbox_to_panel(pannello_Hermite,"Modifica derivate",&modifica_derivata);

	glui ->add_column(FALSE);

	//Bezier
	pannello_Bezier = glui->add_panel("Bezier", GLUI_PANEL_EMBOSSED);
	buttonBezier = glui->add_button_to_panel(pannello_Bezier,"BEZIER",2,scelta_metodi);

	//possibilità di eseguire subdivision
	glui->add_checkbox_to_panel(pannello_Bezier,"Attiva Suddivisione",&attiva_suddivisione);

	//valore del parametro t per la quale valutare la subdivision
	spinner_subd = glui->add_spinner_to_panel(pannello_Bezier, "Suddivisione per ", GLUI_SPINNER_FLOAT, &valore_suddivisione);
	spinner_subd -> set_speed(0.2);

	glui ->add_column(FALSE);

	//Spline
	pannello_Spline = glui->add_panel("Spline", GLUI_PANEL_EMBOSSED);
	bottone_Spline = glui->add_button_to_panel(pannello_Spline,"SPLINE",3,scelta_metodi);

	radio_parametrizzazione = glui->add_radiogroup_to_panel(pannello_Spline, &scelta_parametrizzazione_spline);
	glui->add_radiobutton_to_group(radio_parametrizzazione, "Parametrizzazione Uniforme");
	glui->add_radiobutton_to_group(radio_parametrizzazione, "Parametrizzazione delle Corde");

	//Modifica Molteplicità
	glui->add_checkbox_to_panel(pannello_Spline,"Modifica Molteplicita'",&modifica_molteplicità);
	spinner_i_nodo = glui -> add_spinner_to_panel(pannello_Spline, "Nodo da modificare", GLUI_SPINNER_INT, &indice_nodo);
	spinner_i_nodo -> set_speed(0.1);
	spinner_molteplicità = glui -> add_spinner_to_panel(pannello_Spline, "Molteplicita'", GLUI_SPINNER_INT, &valore_molteplicità);
	spinner_molteplicità -> set_speed(0.1);
	spinner_molteplicità -> set_int_limits(1,ordine_Spline);
	glui->set_main_gfx_window(winIdPrincipale);

	glui ->add_column(FALSE);

	//Pannello Opzioni
	pannello_opzioni = glui->add_panel("Opzioni", GLUI_PANEL_EMBOSSED);
	radio_opzioni = glui->add_radiogroup_to_panel(pannello_opzioni, &scelta_opzioni);
	glui->add_radiobutton_to_group(radio_opzioni, "Aggiunta punti");
	glui->add_radiobutton_to_group(radio_opzioni, "Modifica punti");
	glui->add_radiobutton_to_group(radio_opzioni, "Modifica pesi");
	glui->add_button_to_panel(pannello_opzioni, "Punti Preimpostati",4,scelta_metodi);

}

void reshape(GLsizei width, GLsizei height){
	//metodo per non fare modificare la grandezza della finestra principale
	if(width != larghezza_principale || height != altezza_principale)
	{
		glutReshapeWindow(larghezza_principale, altezza_principale);
	}
}

void main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(larghezza_principale,altezza_principale);
	glutInitWindowPosition(1,1);
	winIdPrincipale = glutCreateWindow("Curve Matematiche 2D");

	//per impostare l'icona nella finestra principale
	HWND hwnd = FindWindow(NULL, _T("Curve Matematiche 2D") );
	HANDLE icon = LoadImage(GetModuleHandle(NULL), _T("Icon.ico"), IMAGE_ICON, 64, 64, LR_LOADFROMFILE | LR_COLOR);
	SendMessage(hwnd, (UINT)WM_SETICON, ICON_BIG, (LPARAM)icon);

	glutDisplayFunc(display);
	glutMouseFunc(myMouse);
	glutMotionFunc(mouseMove);
	glutReshapeFunc(reshape); //per fare in modo che la finstra principale non venga ridimensionata

	glutInitWindowSize(larghezza_secondaria,altezza_secondaria);
	glutInitWindowPosition(30 + larghezza_principale,1);
	winIdFunzioniBase = glutCreateWindow("Funzioni Base");
	glutDisplayFunc(display);

	myinit();

	//per impostare l'icona nella finestra delle funzioni base
	hwnd = FindWindow(NULL, _T("Funzioni Base") );
	SendMessage(hwnd, (UINT)WM_SETICON, ICON_BIG, (LPARAM)icon);

	//per creare il pannello delle opzioni nella windows Principale
	createOptionGlui();

	glutMainLoop();
}