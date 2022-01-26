/*
 * voltageDivider.h
 *
 *  Created on: Apr 6, 2021
 *      Author: Sam S. AME
 */

#ifndef VOLTAGEDIVIDER_H_
#define VOLTAGEDIVIDER_H_
#include <WDF_Creator_core1/src/wdf.h>
#include "common/audio_system_config.h"

//#define PI 3.14159265358

class voltageDivider : wdfTree
{
private:
	Capacitor* C1;
	Resistor* R2;
	Series* S1;
//	Inverter* I1;
	IdealVoltageSource* Vs;

public:
	voltageDivider()
	{
		setSampleRate(AUDIO_SAMPLE_RATE);
		setAlpha(1.0);

		C1 = new Capacitor(1e-6,sr,alpha);
		R2 = new Resistor((float)1000.0);
		S1 = new Series (C1, R2);
//		I1 = new Inverter(S1);
		Vs = new IdealVoltageSource();

		S1->connectToNode(Vs);
	}

	void setInputVoltage(float inSamp)
	{
		Vs->setVoltage(inSamp);
	}

	float getOutputVoltage()
	{
		return C1->getPortVoltage();
	}

	float processSample(float inSamp)
	{
		Vs->setVoltage(inSamp);
		Vs->calcIncidentWave(S1->calcReflectedWave()); //calculate waves incident to root by reflecting from leaves
		float outSamp = C1->getPortVoltage();
		S1->calcIncidentWave(Vs->calcReflectedWave()); //calculate waves reflected from root incident to leaves

		return outSamp;
	}

	void setCutoffFreq(float fc)
	{
		float C = 1e-6;
		float Res = 1.0/(2.0 * PI * C * fc);
		R2->setResistance(Res);
	}


};



#endif /* VOLTAGEDIVIDER_H_ */
