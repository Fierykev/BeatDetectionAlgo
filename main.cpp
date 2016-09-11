#define _WIN32_DCOM
#define _CRT_SECURE_NO_DEPRECATE
#define _USE_MATH_DEFINES

#include <windows.h>
#include <xaudio2.h>
#include <strsafe.h>
#include <shellapi.h>
#include <mmsystem.h>
#include <conio.h>
#include <math.h>
#include <inttypes.h>

#include <limits>

#include <bitset>

#include <vector>

#include "SDKwavefile.h"

#include "FFT.h"

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

// TMP

#include <iostream>
#include <fstream>

using namespace std;

// Microsoft Helper Macros

#ifndef SAFE_DELETE_ARRAY
#define SAFE_DELETE_ARRAY(p) { if(p) { delete[] (p);   (p)=NULL; } }
#endif
#ifndef SAFE_RELEASE
#define SAFE_RELEASE(p)      { if(p) { (p)->Release(); (p)=NULL; } }
#endif

#define swap(x, y)temp = x; x = y; y = temp; // temp is found in the FFT method

#define Swap(x, y)wtemp = x; x = y; y = wtemp; // temp is found in the FFT method

// definitions

//#define SUBBANDS 32
/*GRANDE
#define SONGSCALE .5

#define ZOOM 100.0

#define SMOOTH 1//1000

#define LINEARINCA 0//1.0

#define LINEARINCB 64//358.0 / 36.0

#define CENERGY 7//7.5//15 / 2 

#define CVARIENCE pow(10, 11)

#define SCREENW 1500

#define SCREENH 600
*/

#define SONG L"Audio/Angel_Beats_ED Mono 16.wav"

#define SONGSCALE 1.0

#define ZOOM 100.0

#define SMOOTH 10

#define NUMCHECKS 2

#define LINEARINCA 0//1.0

#define LINEARINCB 32//358.0 / 36.0

double CENERGY = 9;//25/2//7.5//15 / 2

#define CENERGY2 15

#define CVARIENCE 0//pow(10,8)//(pow(10, 8) + 500000)

#define BEATSEARCHTHRESH 9216

#define SCREENW 1500

#define SCREENH 600

// Beat Detection Globals

#define SCALE_E .001

#define beatmin 6

#define beatmax 14

#define MINES 1000

#define BEATSAMPLE 1024.0

double* Ei, *Es, *V;

double subband = 0;

double subbanddatasize = 0;

// Window Draw Globals

double movemulti = 0; // shift screen over

double linepos = 0;

double beatpos = 0;

double pageend = 0;

void UpdateWindow(sf::RenderWindow& App, XAUDIO2_VOICE_STATE& state, vector<int>* beats, int size)
{
	// draw image

	App.clear();

	sf::RectangleShape line(sf::Vector2f(1, 200));

	linepos = ZOOM * (float)state.SamplesPlayed / (float)44100 - SCREENW * movemulti;

	if (SCREENW < linepos)
	{
		// find new pageend

		for (int i = pageend; i < beats->size(); i++)
		{
			beatpos = ZOOM * (((float)beats->at(i)) * size / 44100.0) - SCREENW * movemulti;

			if (SCREENW < beatpos)
			{
				pageend = i;

				break;
			}
		}

		// calculate linepos

		movemulti++;

		linepos = ZOOM * (float)state.SamplesPlayed / (float)44100 - SCREENW * movemulti;
	}

	line.setPosition(linepos, 50);

	App.draw(line);

	// if draw it with line follow

	for (int i = pageend; i < beats->size(); i++)
	{
		sf::RectangleShape shape(sf::Vector2f(1, 100));

		beatpos = ZOOM * (((float)beats->at(i)) * size / 44100.0) - SCREENW * movemulti;

		if (SCREENW < beatpos)
			break;

		if (beatpos < linepos)
			shape.setFillColor(sf::Color(255, 0, 0, 255));
		else
			shape.setFillColor(sf::Color(255, 255, 255, 255));

		shape.setPosition(beatpos, 50);

		App.draw(shape);
	}

	const int nlines = SCREENW / (ZOOM / subband * size / 44100.0);

	// draw graph for Es

	sf::VertexArray linesEs(sf::LinesStrip, nlines);

	for (int i = 0; i < nlines && i + movemulti * nlines < subbanddatasize; i++)
	{
		linesEs[i] = sf::Vector2f(ZOOM * (i / subband * size / 44100.0), SCREENH - SCALE_E * Es[(int)(i + movemulti * nlines)]);

		if (linesEs[i].position.x < linepos)
			linesEs[i].color = sf::Color(255, 255, 0, 255);
		else
			linesEs[i].color = sf::Color(0, 255, 0, 255);
	}

	App.draw(linesEs);

	// draw graph for Ei

	sf::VertexArray linesEi(sf::LinesStrip, nlines);

	for (int i = 42; i < nlines && i + movemulti * nlines < subbanddatasize; i++)
	{
		linesEi[i] = sf::Vector2f(ZOOM * (i / subband * size / 44100.0), SCREENH - SCALE_E * Ei[(int)(i + movemulti * nlines) - 42] * CENERGY);

		if (linesEi[i].position.x < linepos)
			linesEi[i].color = sf::Color(255, 0, 0, 255);
	}

	App.draw(linesEi);
	/*
	// TMP

	sf::RectangleShape lineb(sf::Vector2f(1, 100));

	lineb.setPosition(ZOOM * (1024.0 / subband * size / 44100.0), 50);

	lineb.setFillColor(sf::Color(0, 255, 0, 255));

	App.draw(lineb);
	*/
	/*
	// draw graph for V

	sf::VertexArray linesV(sf::LinesStrip, nlines);

	for (int i = 42; i < nlines && i + movemulti * nlines < subbanddatasize; i++)
	{
	linesV[i] = sf::Vector2f(ZOOM * (i / subband * size / 44100.0), SCREENH - .000001 * V[(int)(i + movemulti * nlines)]);

	if (linesV[i].position.x < linepos)
	linesV[i].color = sf::Color(0, 255, 255, 255);
	else
	linesV[i].color = sf::Color(0, 0, 255, 255);
	}

	App.draw(linesV);
	*/
	// show the window

	App.display();

	/*
	for (int i = 0; i < beats->size(); i++)
	{
	sf::RectangleShape shape(sf::Vector2f(1, 100));

	shape.setPosition(100.0 * ((((float)beats->at(i)) * size / 44100.0 - (float)state.SamplesPlayed / (float)44100)), 50);

	App.draw(shape);
	}*/
}

HRESULT PlayWave(IXAudio2* pXaudio2, CWaveFile& wav, BYTE* pbWaveData, sf::RenderWindow& App, vector<int>* beats, int size)
{
	system("PAUSE");
	///// tmp

	HRESULT hr = S_OK;

	// Get format of wave file
	WAVEFORMATEX* pwfx = wav.GetFormat();

	// Calculate how many bytes and samples are in the wave
	DWORD cbWaveSize = wav.GetSize();

	//
	// Play the wave using a XAudio2SourceVoice
	//

	// Create the source voice
	IXAudio2SourceVoice* pSourceVoice;
	if (FAILED(hr = pXaudio2->CreateSourceVoice(&pSourceVoice, pwfx)))
	{
		wprintf(L"Error %#X creating source voice\n", hr);
		SAFE_DELETE_ARRAY(pbWaveData);
		return hr;
	}

	// Submit the wave sample data using an XAUDIO2_BUFFER structure
	XAUDIO2_BUFFER buffer = { 0 };
	buffer.pAudioData = pbWaveData;
	buffer.Flags = XAUDIO2_END_OF_STREAM;  // tell the source voice not to expect any data after this buffer
	buffer.AudioBytes = cbWaveSize;

	if (FAILED(hr = pSourceVoice->SubmitSourceBuffer(&buffer)))
	{
		wprintf(L"Error %#X submitting source buffer\n", hr);
		pSourceVoice->DestroyVoice();
		SAFE_DELETE_ARRAY(pbWaveData);
		return hr;
	}

	hr = pSourceVoice->Start(0);

	// Let the sound play
	BOOL isRunning = TRUE;
	while (SUCCEEDED(hr) && isRunning)
	{
		XAUDIO2_VOICE_STATE state;
		pSourceVoice->GetState(&state);
		isRunning = (state.BuffersQueued > 0) != 0;

		UpdateWindow(App, state, beats, size);

		// Wait till the escape key is pressed
		if (GetAsyncKeyState(VK_ESCAPE))
			break;

		//Sleep(100);
	}

	pSourceVoice->DestroyVoice();
	SAFE_DELETE_ARRAY(pbWaveData);

	return hr;
}

vector<int>* AnalyzeWave(CWaveFile& wav, BYTE* pbWaveData, int size)
{
	// Based on the http://archive.gamedev.net/archive/reference/programming/features/beatdetection/ article.

	// ONLY FOR: 16 bit mono wave files right now

	// convert the bytes to doubles for the FFT algorithm

	const int samplerate = 32767;

	//const int size = pow(2,9); // pow(2, ceil(log2(wav.GetSize() / 2.0)));

	const int smooth = pow(2, 1);//4096; // smooth the curve by this value

	//const double C = 1.3; // constant for beat detection

	const double bytesize = wav.GetSize(); // length of the byte array from the wave file

	const int datasize = ceil((bytesize / 2.0) / 1024.0) * 1024.0 * SONGSCALE; // the size of the data after converted from bytes

	Complex* datac = new Complex[1024];//datasize];

	Complex* datac2 = new Complex[datasize];

	//bool* energypeak = new bool[(int)bytesize / size + 21];

	std::vector <int>* beats = new vector <int>;

	// load the data to be processed

	for (int l = 0; l < datasize / 1024; l++)
	{
		for (unsigned int i = 0; i < 1024; i++)
		{
			if ((l * 1024 + i) * 2 < bytesize)
				datac[i].real((__int16)((__int16)((pbWaveData[(l * 1024 + i) * 2 + 1] << 8) | pbWaveData[(l * 1024 + i) * 2]) / SMOOTH));
			else
				datac[i].real(0);

			datac[i].imag(0);
		}

		CArray data2(datac, 1024);

		fft(data2);

		for (int i = 0; i < 1024; i++)
			datac2[l * 1024 + i] = data2[i];
	}

	double* B = new double[datasize]; // create double to store the amps

	for (unsigned long i = 0; i < datasize; i++)
	{
		B[i] = sqrt((datac2[i].real() * datac2[i].real()) + (datac2[i].imag() * datac2[i].imag())); // mod(x + iy) = sqrt(x^2 + y^2)
	}

	// calculate the number of subbands

	double sum = 1024;

	int SUBBANDS;

	for (SUBBANDS = 1; 0 < sum; SUBBANDS++)
		sum -= (int)(LINEARINCA * SUBBANDS + LINEARINCB);
	cout << SUBBANDS << ", " << sum << endl;

	SUBBANDS -= 1;

	subband = SUBBANDS;

	Es = new double[datasize / 1024 * SUBBANDS];

	int wi; // width of subband

	// analize for 1024 samples

	for (int l = 0; l < datasize / 1024; l++)
	{
		sum = 0;

		for (int i = 0; i < SUBBANDS; i++) // divide into subbands
		{
			wi = i < SUBBANDS - 1 ? LINEARINCA * (i + 1) + LINEARINCB : 1024 - sum; // width of the subband

			Es[l * SUBBANDS + i] = 0; // clear the value

			for (int k = sum; k < sum + wi; k++)
			{
				Es[l * SUBBANDS + i] += B[l * 1024 + k];
			}

			sum += wi; // add the width to the sum

			Es[l * SUBBANDS + i] *= wi / 1024.0;
		}
	}

	// average subband energy

	subbanddatasize = datasize / 1024 * SUBBANDS;

	Ei = new double[datasize / 1024 * SUBBANDS];

	V = new double[datasize / 1024 * SUBBANDS];

	sum = 0;

	for (int i = 0; i < 43; i++)
		sum += Es[i];

	Ei[0] = sum / 43.0;

	for (int i = 1; i < datasize / 1024 * SUBBANDS - 43; i++)
	{
		// energy calc

		sum += Es[i + 42] - Es[i - 1];

		Ei[i] = sum / 43.0;
	}
	/*
	// varience

	double sumV = 0;

	for (int i = 43; i < datasize / 1024 * SUBBANDS - 43; i+=NUMCHECKS)
	{
	// varience calc

	sumV = 0;

	for (int k = 0; k < 43; k++)
	sumV += (Es[i - 42 + k] - Ei[i]) * (Es[i - 42 + k] - Ei[i]);

	V[i] = sumV / 43.0;
	}*/

	// * 16 angel beats
	// 25 at 54 subbands

	//CENERGY = Es[max + 42] / Ei[max] * 10.0;
	/*
	double* C = new double[datasize / 64];

	int max;

	int CSCALE = 1024;

	for (int l = 0; l < datasize / 1024 * SUBBANDS; l += CSCALE)
	{
	max = l;

	for (int i = 0; i < CSCALE; i++)
	if (Ei[max] < Ei[i + l])
	max = i + l;

	if (max + 42 < datasize / 1024 * SUBBANDS)
	C[l / CSCALE] = Es[max + 42] / Ei[max] * 15.0;
	else
	C[l / CSCALE] = 7.0;
	}

	//printf("%f", (float)CENERGY);

	sum = 0;

	for (int i = 42; i < datasize / 1024 * SUBBANDS; i += NUMCHECKS)
	if (Ei[i - 42] * C[(i - 42) / CSCALE] < Es[i])// && V[i] > CVARIENCE)
	{
	// add the beats to the list

	/*if (BEATSEARCHTHRESH < sum - i)
	{*/
	//beats->push_back(i / SUBBANDS);

	//sum = i;
	/*}
	else // go back and reevaluate the beats and a lower threshold
	{
	for (int j = sum; j < i; j++)
	{
	if (Ei[j - 42] * CENERGY2 < Es[j])// && V[i] > CVARIENCE)
	beats->push_back(j / SUBBANDS);
	}

	beats->push_back(i / SUBBANDS);

	sum = i;
	}*//*
	}*/

	double max = 0;// , maxEs = 0;

	for (int l = 42; l + BEATSAMPLE < datasize / 1024 * SUBBANDS; l += BEATSAMPLE)
	{
		max = 0;

		//maxEs = 0;

		for (int j = 0; j < BEATSAMPLE; j += NUMCHECKS)
		{
			if (max < Es[l + j] / Ei[l + j - 42])
				max = Es[l + j] / Ei[l + j - 42];

			/*if (maxEs < Es[l + j])
			maxEs = Es[l + j];*/
		}

		// round to .5

		max = floor(max);//(int)(max * 10.0) % 10 < 5 ? floor(max) : floor(max) + .5;

		// constrain the bands

		if (max < beatmin)
			max = beatmin;
		else if (beatmax < max)
			max = beatmax;

		for (int j = 0; j < BEATSAMPLE; j += NUMCHECKS)
		{
			if (Ei[l + j - 42] * max < Es[l + j])
				beats->push_back((l + j) / SUBBANDS);
		}

		//cout << (float)maxEs << endl;
	}
	/*
	// filter through the beats and see if there are patterns that need to be fixed

	double deltabeats;

	double pastbeat = 0;

	double deltabeatsh = 0;

	for (int i = 1; i < beats->size(); i++)
	{
	deltabeats = beats->at(i) - pastbeat;

	pastbeat = beats->at(i);

	if (abs(deltabeatsh - deltabeats) < 5)
	{
	beats->erase(beats->begin() + i);

	i--;
	}
	else
	deltabeatsh = deltabeats;
	}
	*/
	// clear heap

	//free(Ei);

	//free(V);

	//free(Es);

	//free(bandsi);

	free(datac);

	free(datac2);

	free(B);

	return beats;
}

HRESULT LoadWave(IXAudio2*& pXaudio2, LPCWSTR szFilename, CWaveFile& wav, BYTE*& pbWaveData)
{
	HRESULT hr = S_OK;

	// locate the wave file

	WCHAR strFilePath[MAX_PATH];

	wcscpy_s(strFilePath, MAX_PATH, szFilename);

	// TODO: ADD A FILE CHECKER
	/*if (FAILED(hr = FindMediaFileCch(strFilePath, MAX_PATH, szFilename)))
	{
	wprintf(L"Failed to find media file: %s\n", szFilename);
	return hr;
	}*/

	// read the wave file

	if (FAILED(hr = wav.Open(strFilePath, NULL, WAVEFILE_READ)))
	{
		wprintf(L"Failed reading WAV file: %#X (%s)\n", hr, strFilePath);
		return hr;
	}

	// Get format of wave file
	WAVEFORMATEX* pwfx = wav.GetFormat();

	// Calculate how many bytes and samples are in the wave
	DWORD cbWaveSize = wav.GetSize();

	// Read the sample data into memory
	pbWaveData = new BYTE[cbWaveSize];

	if (FAILED(hr = wav.Read(pbWaveData, cbWaveSize, &cbWaveSize)))
	{
		wprintf(L"Failed to read WAV data: %#X\n", hr);
		SAFE_DELETE_ARRAY(pbWaveData);
		return hr;
	}

	return hr;
}

HRESULT InitXAudio2(IXAudio2*& pXAudio2)
{
	HRESULT hr;

	// init XAudio2

	CoInitializeEx(NULL, COINIT_MULTITHREADED);

	UINT32 flags = 0;

#ifdef _DEBUG
	flags |= XAUDIO2_DEBUG_ENGINE;
#endif

	// could not create the audio engine

	if (FAILED(hr = XAudio2Create(&pXAudio2, flags)))
	{
		wprintf(L"Failed to init XAudio2 engine: %#X\n", hr);
		CoUninitialize();
		return 0;
	}

	// create a mastering voice

	IXAudio2MasteringVoice* pMasteringVoice = NULL;

	if (FAILED(hr = pXAudio2->CreateMasteringVoice(&pMasteringVoice)))
	{
		wprintf(L"Failed creating mastering voice: %#X\n", hr);
		SAFE_RELEASE(pXAudio2);
		CoUninitialize();
		return 0;
	}
}

int main()
{
	sf::RenderWindow App(sf::VideoMode(SCREENW, SCREENH), "Beat detector");

	IXAudio2* pXAudio2 = NULL;

	InitXAudio2(pXAudio2);

	// finished creating XAudio

	// load a wave

	CWaveFile wav;

	BYTE* pbWaveData;

	LoadWave(pXAudio2, SONG, wav, pbWaveData);

	int size = 1024;

	vector<int>* beats = AnalyzeWave(wav, pbWaveData, size);

	PlayWave(pXAudio2, wav, pbWaveData, App, beats, size);

	free(beats);

	// stop XAudio2

	SAFE_RELEASE(pXAudio2);
	CoUninitialize();

	return 0;
}