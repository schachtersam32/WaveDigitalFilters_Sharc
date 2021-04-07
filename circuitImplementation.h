/*
 * circuitImplementation.h
 *
 *  Created on: Apr 6, 2021
 *      Author: Sam S. AME
 */

#ifndef CIRCUITIMPLEMENTATION_H_
#define CIRCUITIMPLEMENTATION_H_

class LC_Oscillator : wdfTree
{
private:
	Inductor* L;
	Capacitor* C;
	Switch* SW;
	Series* S1;
	Series* S2;
	ResistiveVoltageSource* Vs;

public:
	LC_Oscillator()
	{
		setAlpha(1.0);
		setSampleRate(AUDIO_SAMPLE_RATE);

		C = new Capacitor(1.0e-3,sr,alpha);
		L = new Inductor(1.0e3,sr,alpha);
		Vs = new ResistiveVoltageSource(1.0);

		SW = new Switch;
		S2 = new Series(L,C);
		S1 = new Series(Vs,S2);
		S1->connectToNode(SW);

	}

	~LC_Oscillator(){}
	void setParams(float freq, bool switchClosed)
	{
		float C_val = 1.0e-6;
		float L_val = 1.0/((2*PI*freq)*(2*PI*freq)*C_val);

		C->setCapacitance(C_val);
		L->setInductance(L_val);
		SW->closeSwitch(switchClosed);
	}

	float processSample(float inSamp)
	{
		Vs->setVoltage(inSamp);

		SW->calcIncidentWave(S1->calcReflectedWave());
		S1->calcIncidentWave(SW->calcReflectedWave());
		return C->getPortVoltage();
	}
};



#endif /* CIRCUITIMPLEMENTATION_H_ */
