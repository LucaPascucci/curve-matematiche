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

#pragma comment(linker,"/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

int larghezzaPrincipale = 800, altezzaPrincipale = 800, PuntoSelezionato = -1;
int larghezzaSecondaria = 400, altezzaSecondaria = 300;

//definiamo la nostra interfaccia
int winIdPrincipale, winIdFunzioniBase, scelta_opzioni = 0, mod_molt = 0, i_nodo, molte;
int m = 4; //Ordine spline
GLUI_Panel *pannello_opzioni;
GLUI_RadioGroup *radio_opzioni;

//pannello scelta parametri
int scelta_param = 0;
GLUI_Panel *pannello_param;
GLUI_RadioGroup *radio_param;

//pannello per esecuzione Hermite
GLUI_Panel *pannello_Hermite;
GLUI_Button *buttonHermite;

//pannello per esecuzione Bezier
GLUI_Panel *pannello_Bezier;
GLUI_Button *buttonBezier;

//pannello per esecuzione Spline
GLUI_Panel *pannello_Spline;
GLUI_Button *buttonSpline;

float t_subd;
int alg_subd = 0; //1: viene applicato l'algoritmo subdivision

GLUI_Spinner *spinner_subd, *spinner_i_nodo, *spinner_molte;

//scelta dei metodi
int metodo = 0; //1: Hermite; 2: Bezier; 3: curve Spline

//modifica derivate
int mod_der = 0;

using namespace std;

typedef struct glPoint2D{
	GLfloat x,y;
}GLPOINT2D;

//contenitore di dati
vector <GLPOINT2D> Punti;

//conserva le coordinate dei punti dove è stata modificata la derivata
vector <GLPOINT2D> DerivateMod;

void myMouse(int button, int state, GLint xmouse, GLint ymouse){
	
	int TOLL = 3;
	float distanza, distanza1;

	GLPOINT2D newPoint, zero;
	zero.x = 0.0;
	zero.y = 0.0;
	newPoint.x = xmouse;
	newPoint.y = altezzaPrincipale - ymouse;

	if (state == GLUT_DOWN){ //se lo stato del bottone è premuto
		switch(button){

		case GLUT_LEFT_BUTTON:

			//scelta_opzioni = 0 --> modalita' inserimento
			//scelta_opzioni = 1 --> modifica punti inseriti
			//mod_der = 0 --> se non voglio modificare la derivata sui punti inseriti
			//mod_der = 1 --> se voglio modificare la derivata sui punti inseriti
			if(scelta_opzioni == 1 || mod_der == 1){

				//ho già inserito dei punti, devo "catturare" il punto cliccato
				//tramite un controllo "ad area" che confronta la zona cliccata
				//coi punti inseriti, se la distanza è minore di una costante
				//allora selezionerò il punto per poi modificarne la posizione
				if(Punti.size() > 0){
					PuntoSelezionato = 0;

					//calcoliamo la distanza dal newpoint da tutti i punti inseriti prima
					//e prendo l'indice del punto più vicino
					distanza = sqrt((Punti.at(0).x - newPoint.x)*(Punti.at(0).x - newPoint.x) + (Punti.at(0).y - newPoint.y)*(Punti.at(0).y - newPoint.y));

					for (int i=1; i<Punti.size(); i++){
						distanza1 = sqrt((Punti.at(i).x - newPoint.x)*(Punti.at(i).x - newPoint.x) + (Punti.at(i).y - newPoint.y)*(Punti.at(i).y - newPoint.y));

						//faccio il controllo sui minimi delle distanze
						if (distanza1<distanza){
							PuntoSelezionato = i;
							distanza = distanza1;
						}
					}

					//confronto ora col mio indice di tolleranza
					if(distanza > TOLL){
						PuntoSelezionato = -1;
					}
				}
			} else if (scelta_opzioni == 0){
				glBegin(GL_POINTS);
				glVertex2f(newPoint.x, newPoint.y);
				glEnd();

				glFlush();

				Punti.push_back(newPoint);

				DerivateMod.push_back(zero);
			}
			break;

		case GLUT_RIGHT_BUTTON:
			Punti.clear();
			break;

		case GLUT_MIDDLE_BUTTON:
			Punti.pop_back();
			break;
		}
	} else {
		switch(button){
		case GLUT_LEFT_BUTTON:
			PuntoSelezionato = -1;
			break;
		}

	}

	//forza il ridisegno
	glutPostRedisplay();
}

void mouseMove(GLint xmouse, GLint ymouse){

	GLPOINT2D newPoint;
	newPoint.x = xmouse;
	newPoint.y = altezzaPrincipale - ymouse;

	if(PuntoSelezionato >= 0){
		if(scelta_opzioni == 1){
			Punti.at(PuntoSelezionato) = newPoint;
		} else if (mod_der == 1){
			//la moltiplicazione per 5 è una regola
			float derx = (newPoint.x - Punti.at(PuntoSelezionato).x)*5;
			float dery = (newPoint.y - Punti.at(PuntoSelezionato).y)*5;

			DerivateMod.at(PuntoSelezionato).x = derx;
			DerivateMod.at(PuntoSelezionato).y = dery;
		}
	}

	glutPostRedisplay();
}

void parametrizzazione_uniforme(float* t){

	//definisco il passo della parametrizzazione, dividendo l'ampiezza dell'intervallo
	//per il numero di punti presenti (l'ampiezza dell'intervallo è sempre 1)
	float step = 1.0/(float)(Punti.size()-1);

	for (int i = 0; i<Punti.size(); i++){
		t[i] = i*step;
	}

}

void parametrizzazione_corde(float* t){

	t[0] = 0;

	for (int i = 1; i<Punti.size(); i++){
		t[i] = t[i-1] + sqrt((Punti.at(i).x - Punti.at(i-1).x)*(Punti.at(i).x - Punti.at(i-1).x) + (Punti.at(i).y - Punti.at(i-1).y)*(Punti.at(i).y - Punti.at(i-1).y));
	}

	//divido per il massimo dell'ultimo componente per riportarlo nel range 0-1
	for (int i = 0; i<Punti.size(); i++){
		t[i] = t[i]/t[Punti.size()-1];
	}

}

//definizione funzioni base di Hermite
#define PHI0(t) (2.0*t*t*t - 3.0*t*t + 1)
#define PHI1(t) (t*t*t - 2.0*t*t + t)

#define PSI0(t) (-2.0*t*t*t + 3.0*t*t)
#define PSI1(t) (t*t*t - t*t)

//derivate rispetto a t della componente parametrica in x

//metodo del rapporto incrementale
float dx(int i, float* t){
	if (i <= 0)
		return 0;

	//altrimenti restituisco la derivata in x
	return (Punti.at(i).x - Punti.at(i-1).x)/(t[i] - t[i-1]);
}

//derivate rispetto a t della componente parametrica in y

//metodo del rapporto incrementale
float dy(int i, float* t){
	if (i <= 0)
		return 0;

	//altrimenti restituisco la derivata in x
	return (Punti.at(i).y - Punti.at(i-1).y)/(t[i] - t[i-1]);
}

float DX(int i, float *t){
	
	if(DerivateMod.at(i).x == 0)
		return dx(i, t);
	
	if(DerivateMod.at(i).x != 0)
		return DerivateMod.at(i).x;
}

float DY(int i, float *t){
	
	if(DerivateMod.at(i).y == 0)
		return dy(i, t);
	
	if(DerivateMod.at(i).y != 0)
		return DerivateMod.at(i).y;
}

//implementazione interpolazione di Hermite
void InterpolazioneHermite(float* t){

	//valutiamo la nostra curva su 1000 (mila) valori
	//numero di valori del parametro t in cui valutare la curva interpolante di Hermite
	int nvpt = 1000;

	float passot = 1.0/(float)(nvpt-1);

	//per ogni valore, devo vedere in quale sottointervallo cade il mio punto di valutazione
	int is = 0; //Indice del sottointervallo a cui appartiene il valore del parametro in cui valutare la curva
	
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0,0.4,0.0);
	for (float ti = 0; ti <= 1; ti += passot){

		if (ti > t[is+1])
			is++; //is è l'indice dell'intervallo a cui appartiene il valore ti del parametro

		//mappare il valore ti (appartenente all'intervallo [t[is],t[is+1]) nell'intervallo [0,1]
		//perchè le funzioni base sono state definite su un intervallo [0,1]
		float tim;
		float range = (t[is+1] - t[is]);
		tim = (ti - t[is])/range;

		GLPOINT2D IH;

		IH.x = Punti.at(is).x*PHI0(tim) + DX(is,t)*PHI1(tim)*range + Punti.at(is+1).x*PSI0(tim) + DX(is+1,t)*PSI1(tim)*range;
		IH.y = Punti.at(is).y*PHI0(tim) + DY(is,t)*PHI1(tim)*range + Punti.at(is+1).y*PSI0(tim) + DY(is+1,t)*PSI1(tim)*range;

		glVertex2f(IH.x, IH.y);
	}
	glEnd();
}

//algoritmo di subdivision
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
			c[i].x = c[i].x * (1 - t_subd) + t_subd * c[i+1].x;
			c[i].y = c[i].y * (1 - t_subd) + t_subd * c[i+1].y;
		}

		c1[j].x = c[0].x;
		c1[j].y = c[0].y;
		c2[j] = c[Punti.size() - 1 - j];
		//c2[j].y = c[Punti.size() - 1 - j].y;
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

	int npv = 1000;
	float tstep = 1.0/(float)(npv - 1);

	GLPOINT2D *c = new GLPOINT2D[Punti.size()];

	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 0.0, 0.0);
	for(float ti = 0; ti <= 1; ti += tstep){
		for(int i = 0; i < Punti.size() ; i++)
		{
			c[i] = Punti.at(i);
			//c[i] = Pesi[i] * Punti.at(i);
			//w[i] = Pesi[i];
		}

		for(int j = 1; j <= Punti.size(); j++){
			for(int i = 0; i < Punti.size() - j ; i++){
				c[i].x = c[i].x * (1 - ti) + ti * c[i+1].x;
				c[i].y = c[i].y * (1 - ti) + ti * c[i+1].y;
				//Pesi[i] = Pesi[i].y * (1 - ti) + ti * Pesi[i+1].y;
				//w[i] = w[i] * (1 - ti) + ti * w[i+1];
			}
		}
		glVertex2f(c[0].x, c[0].y);
		//glVertex2f(c[0].x / Pesi[0], c[0].y / Pesi[0]);
		//glVertex2f(c[0].x / w[0], c[0].y / w[0]);
	}
	glEnd();

	delete(c);
}

void disegnaBezier(){

	int campioni = 200;
	float ti = 0, dt = 1.0/(campioni - 1);

	//allochiamo la matrice B, nella quale salveremo le nostre funzioni base
	float **B = new float *[campioni];

	//posizione dove memorizzare la valutazione che corrisponde al valore del parametro t
	for (int k = 0; k < campioni; k++, ti += dt){
		B[k] = new float[Punti.size()+1]();
		B[k][Punti.size()] = 1.0; //condizione iniziale

		//formula ottimizzata, per evitare controlli su valori nulli
		//ricorsiva, ci salviamo solo il valore che ci interessa
		for (int i = 1; i < Punti.size(); i++){
			float d1b = 0;
			for (int j = 0; j < Punti.size(); j++){
				B[k][j] = (1-ti)*B[k][j+1]+d1b;
				d1b = ti*B[k][j+1];
			}
			B[k][Punti.size()] = d1b;
		}
	}

	//disegnamo le funzioni base, che sono le colonne di questa matrice B
	for (int ibase = 0; ibase < Punti.size(); ibase++){
		glBegin(GL_LINE_STRIP);
		for (int l = 0; l < campioni; l++){
			glVertex2f(dt*(float)l, B[l][ibase+1]);
		}
		glEnd();
	}
}

void costruisci_Nodi(float *t, float *Nodi, char* molt)
{
	int i, cont;
	int k = Punti.size() - m; //Numero di Nodi interni all'intervallo

	//Nodi fittizi a sinistra
	for (i = 0; i < m; i++)
	{
		Nodi[i] = 0;
		molt[i] = '4';
	}

	//Costruzione nodi veri
	cont = 2;
	for (i = m; i < m + k; i++)
	{
		Nodi[i] = t[cont];
		molt[i] = '1';
		cont++;
	}

	//Nodi fittizi a destra
	for (i = m + k; i < 2 * m + k; i++)
	{
		Nodi[i] = 1;
		molt[i] = '4';
	}

	if (mod_molt == 1)
	{
		float val_nodo = Nodi[i_nodo];
		if (molte == 2)
		{
			molt[i_nodo] = '2';
		}
		if (molte == 3)
		{
			molt[i_nodo] = '3';
		}
		if (molte > 1)
		{
			for (i = 1; i < molte; i++)
			{
				Nodi[i_nodo + i -1] = val_nodo;
			}
		}
	}
}

int localizza_intervallo_internodale(float t, float *Nodi)
{
	//Implementazione del metoo di bisezione
	int a = m - 1;
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

void disegnaBaseSpline(float *Nodi, char *molt)
{
	int Ncampioni = 200;
	float ti = 0;
	float dt = 1.0/(float)(Ncampioni - 1);

	float **B = new float*[Punti.size()];
	for (int i = 0; i < Punti.size(); i++) {
		B[i] = new float[Ncampioni]();
	}
	//valutazione di ciascuna funzione base per ogni valore del parametro t mediante le formule di Cocks
	for (int k = 0; k < Ncampioni; k++, ti += dt)
	{
		int l = localizza_intervallo_internodale(ti, Nodi);

		B[l][k] = 1;
		for (int i = 0; i < m - 1; i++)
		{
			float tmp = 0.0;
			for (int j = l - i; j <= l; j++)
			{
				float d1 = ti -Nodi[j];
				float d2 = Nodi[i+j+1] - ti;
				float beta = B[j][k] / (d1 + d2);
				B[j-1][k] = d2 * beta + tmp;
				tmp = d1 * beta;
			}
			B[l][k] = tmp;
		}
	}
	for (int i = 0; i < Punti.size(); i++)
	{
		ti = 0; 
		glBegin(GL_LINE_STRIP);
		for (int k = 0; k < Ncampioni; k++, ti += dt)
		{
			glVertex2f(ti, B[i][k]);
		}
		glEnd();
	}
	for (int j = 0; j < Punti.size() + m; j++)
	{
		glColor3f(0.0, 1.0, 0.0);
		glRasterPos2f(Nodi[j], 0.0);
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, molt[j]);
	}
}

void DeBoor(float *t, float *Nodi)
{
	int nvalorit = 1000, i;
	float tstep = 1.0 / (float)(nvalorit - 1);

	GLPOINT2D *c = new GLPOINT2D[m];

	glBegin(GL_LINE_STRIP);
	for (float vt = 0; vt <= 1; vt += tstep)
	{
		int l = localizza_intervallo_internodale(vt, Nodi);
		//Implementamo l'algoritmo di DeBoor
		for (i = 0; i < m; i++)
		{
			c[i] = Punti.at(i + l - m + 1);
		}
		for (int j = 0; j < m - 1; j++)
		{
			for (i = m - 1; i > j; i--)
			{
				int ti = l - m + 1 + i;
				int timj = l + i - j;

				float den = Nodi[timj] - Nodi[ti];
				float dt = (vt - Nodi[ti]) / den;
				
				c[i].x = c[i].x * dt + (c[i-1].x * (1 - dt));
				c[i].y = c[i].y * dt + (c[i-1].y * (1 - dt));
			}
		}
		glVertex2f(c[m - 1].x, c[m - 1].y);
	}
	glEnd();
}

void scelta_metodi(int scelta){

	switch(scelta){
	case 1:
		metodo = 1;
		glutPostRedisplay();
		break;
	
	case 2:
		metodo = 2;
		glutPostRedisplay();
		break;

	case 3:
		metodo = 3;
		glutPostRedisplay();
		break;
	}

}

void display(){

	glutSetWindow(winIdFunzioniBase);
	glClear(GL_COLOR_BUFFER_BIT); //pulisco la finestra secondaria (delle funzioni base)

	glutSetWindow(winIdPrincipale);
	glClear(GL_COLOR_BUFFER_BIT);
	//Definisco il sistema di riferimento per la finestra principale di disegno delle curve
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,float(larghezzaPrincipale),0.0,float(altezzaPrincipale));

	glColor3f(0.0, 0.0, 0.0); //disegno i punti della
	glBegin(GL_POINTS);
	for (int i = 0; i < Punti.size(); i++)
		glVertex2f(Punti.at(i).x, Punti.at(i).y);
	glEnd();
	
	glColor3f(1.0,0.0,0.0);
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < Punti.size(); i++)
		glVertex2f(Punti.at(i).x, Punti.at(i).y);
	glEnd();
	glFlush();

	//Parametrizzazione

	float* t = new float[Punti.size()];

	if (scelta_param == 0){
		parametrizzazione_uniforme(t);
	} else if (scelta_param == 1){
		parametrizzazione_corde(t);
	}

	if(Punti.size() > 1){

		//curve interpolanti di Hermite
		if (metodo == 1){
			InterpolazioneHermite(t);

			if(PuntoSelezionato >= 0 && mod_der == 1){
				//il punto selezionato
				GLPOINT2D P0 = Punti.at(PuntoSelezionato);

				//coordinate del punto corrispondente al valore del parametro t = 0
				//appartenente alla retta tangente alla curva nel punto selezionato
				GLPOINT2D P1;

				P1.x = P0.x - DerivateMod.at(PuntoSelezionato).x * t[PuntoSelezionato];
				P1.y = P0.y - DerivateMod.at(PuntoSelezionato).y * t[PuntoSelezionato];
			
				//coordinate del punto corrispondente al valore del parametro t = 1
				//appartenente alla retta tangente alla curva nel punto selezionato
				GLPOINT2D P2;

				P2.x = P0.x + DerivateMod.at(PuntoSelezionato).x - DerivateMod.at(PuntoSelezionato).x * t[PuntoSelezionato];
				P2.y = P0.y + DerivateMod.at(PuntoSelezionato).y - DerivateMod.at(PuntoSelezionato).y * t[PuntoSelezionato];

				glBegin(GL_LINE_STRIP);
					glVertex2f(P1.x, P1.y);
					glVertex2f(P2.x, P2.y);
				glEnd();
			}
		} else if (metodo == 2){
			
			Bezier();
			if(alg_subd == 1){
				Subdivision();
			}
			glutSetWindow(winIdFunzioniBase);
			glColor3f(1.0, 0.0, 0.0);
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();

			disegnaBezier();
			glMatrixMode(GL_MODELVIEW);

		} else if (metodo == 3) {
			if (Punti.size() >= m) {
				float *Nodi = new float[Punti.size() + 2 * m];
				char *molt = new char[Punti.size() + m];
				costruisci_Nodi(t, Nodi, molt);
				DeBoor(t, Nodi);

				glutSetWindow(winIdFunzioniBase);
				glColor3f(1.0, 0.0, 0.0);
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				disegnaBaseSpline(Nodi, molt);
			}
		}
	}

	glFlush();
}

void myinit (void)
{
	glutSetWindow(winIdFunzioniBase);
	glClearColor(1.0, 1.0, 1.0, 0.0); //colore dello sfondo finestra secondaria
	gluOrtho2D(-0.05, 1.05, -0.05, 1.05);
	glutSetWindow(winIdPrincipale);
	glClearColor(1.0, 1.0, 1.0, 0.0); //colore dello sfondo finestra principale
	glPointSize(6.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,float(larghezzaPrincipale),0.0,float(altezzaPrincipale));
}

void createGlui(){

	GLUI *glui = GLUI_Master.create_glui("Opzioni", GLUI_SUBWINDOW_TOP, 30 + larghezzaPrincipale, 80 + altezzaSecondaria);
	pannello_opzioni = glui->add_panel("Opzioni", GLUI_PANEL_EMBOSSED);
	radio_opzioni = glui->add_radiogroup_to_panel(pannello_opzioni, &scelta_opzioni);
	glui->add_radiobutton_to_group(radio_opzioni, "Aggiunta punti");
	glui->add_radiobutton_to_group(radio_opzioni, "Modifica punti");

	pannello_param = glui->add_panel("Parametrizzazione", GLUI_PANEL_EMBOSSED);
	radio_param = glui->add_radiogroup_to_panel(pannello_param, &scelta_param);
	glui->add_radiobutton_to_group(radio_param, "Uniforme");
	glui->add_radiobutton_to_group(radio_param, "Corde");

	pannello_Hermite = glui->add_panel("Hermite", GLUI_PANEL_EMBOSSED);
	//quando viene premuto il pulsante, viene eseguita l'opzione 1 del metodo scelta_metodi
	//1 è il valore che viene passato alla funzione scelta_metodi
	buttonHermite = glui->add_button_to_panel(pannello_Hermite,"HERMITE",1,scelta_metodi);

	//possibilità di modificare le derivate per modificare il comportamento delle nostre curve interpolate
	glui->add_checkbox_to_panel(pannello_Hermite,"Modifica derivate",&mod_der);

	pannello_Bezier = glui->add_panel("Bezier", GLUI_PANEL_EMBOSSED);
	//quando viene premuto il pulsante, viene eseguita l'opzione 1 del metodo scelta_metodi
	//1 è il valore che viene passato alla funzione scelta_metodi
	buttonBezier = glui->add_button_to_panel(pannello_Bezier,"BEZIER",2,scelta_metodi);

	//possibilità di eseguire subdivision
	glui->add_checkbox_to_panel(pannello_Bezier,"Subdivision",&alg_subd);

	//valore del parametro t per la quale valutare la subdivision
	spinner_subd = glui->add_spinner_to_panel(pannello_Bezier, "Parametro t per Subdivision", GLUI_SPINNER_FLOAT, &t_subd);
	spinner_subd -> set_speed(0.3);

	pannello_Spline = glui->add_panel("Spline", GLUI_PANEL_EMBOSSED);
	buttonSpline = glui->add_button_to_panel(pannello_Spline,"SPLINE",3,scelta_metodi);
	glui->add_checkbox_to_panel(pannello_Spline,"Modifica Molt",&mod_molt);
	spinner_i_nodo = glui -> add_spinner_to_panel(pannello_Spline, "Nodo da modificare", GLUI_SPINNER_INT, &i_nodo);
	spinner_i_nodo -> set_speed(0.1);
	spinner_molte = glui -> add_spinner_to_panel(pannello_Spline, "Molteplicita'", GLUI_SPINNER_INT, &molte);
	spinner_molte -> set_speed(0.1);

	glui->set_main_gfx_window(winIdPrincipale);
}

void reshape(GLsizei width, GLsizei height){
	//metodo per non fare modificare la grandezza della finestra principale
	if(width != larghezzaPrincipale || height != altezzaPrincipale)
	{
		glutReshapeWindow(larghezzaPrincipale, altezzaPrincipale);
	}
}

void main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(larghezzaPrincipale,altezzaPrincipale);
	glutInitWindowPosition(1,1);
	winIdPrincipale = glutCreateWindow("Curve Matematiche 2D");

	//per impostare l'icona nella finestra principale
	HWND hwnd = FindWindow(NULL, _T("Curve Matematiche 2D") );
	HANDLE icon = LoadImage(GetModuleHandle(NULL), _T("Icon.ico"), IMAGE_ICON, 64, 64, LR_LOADFROMFILE | LR_COLOR);
	SendMessage(hwnd, (UINT)WM_SETICON, ICON_BIG, (LPARAM)icon);
	
	glutDisplayFunc(display);
	glutMouseFunc(myMouse);
	glutMotionFunc(mouseMove);
	glutReshapeFunc(reshape); //per fare in modo che la finstra principale non venga resizeta

	glutInitWindowSize(larghezzaSecondaria,altezzaSecondaria);
	glutInitWindowPosition(30 + larghezzaPrincipale,1);
	winIdFunzioniBase = glutCreateWindow("Funzioni Base");
	glutDisplayFunc(display);

	myinit();

	//per impostare l'icona nella finestra delle funzioni base
	hwnd = FindWindow(NULL, _T("Funzioni Base") );
	SendMessage(hwnd, (UINT)WM_SETICON, ICON_BIG, (LPARAM)icon);

	createGlui();
	//per impostare l'icona nella finestra delle opzioni
	hwnd = FindWindow(NULL, _T("Opzioni") );
	SendMessage(hwnd, (UINT)WM_SETICON, ICON_BIG, (LPARAM)icon);


	

	glutMainLoop();
}