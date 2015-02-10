/*
 * Original work Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file VideoNLBayes.cpp
 * @brief Video NL-Bayes denoising functions
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 * @author Pablo Arias <pariasm@gmail.com>
 **/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>

#include "VideoNLBayes.h"
#include "LibMatrix.h"
#include "../Utilities/Utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace VideoNLB
{

/**
 * @brief Initialize Parameters of the NL-Bayes algorithm.
 *
 * @param o_params   : will contain the nlbParams for the first step of the algorithm;
 * @param p_step     : select first or second step;
 * @param p_sigma    : standard deviation of the noise;
 * @param p_size     : size of the video;
 * @param p_flatArea : if true, use the homogeneous area trick for the first step;
 * @param p_verbose  : if true, print some informations.
 * @param p_timeSearchRagneFwd : temporal search range forwards.
 * @param p_timeSearchRagneBwd : temporal search range backwards.
 *
 * @return none.
 **/
void initializeNlbParameters(
	nlbParams &o_params
,	const unsigned p_step
,	const float p_sigma
,	const VideoSize &p_size
,	const bool p_flatArea
,	const bool p_verbose
,	const unsigned timeSearchRangeFwd
,	const unsigned timeSearchRangeBwd
){
	const bool s1 = (p_step == 1);

	//! Standard deviation of the noise
	o_params.sigma = p_sigma;

	//! Size of patches
	if (p_size.channels == 1)
	{
		if(s1) o_params.sizePatch = (p_sigma < 30.f ? 5 : 7);
		else   o_params.sizePatch = 5;
	}
	else
	{
		if(s1) o_params.sizePatch = (p_sigma < 20.f ? 3 : (p_sigma < 50.f ? 5 : 7));
		else   o_params.sizePatch = (p_sigma < 50.f ? 3 : (p_sigma < 70.f ? 5 : 7));
	}

	//! Number of similar patches
	if (p_size.channels == 1)
	{
		if(s1) o_params.nSimilarPatches = (p_sigma < 10.f ? 35 : (p_sigma < 30.f ? 45 : (p_sigma < 80.f ? 90 : 100))); // ASK MARC differs from manuscript
		else   o_params.nSimilarPatches = (p_sigma < 20.f ? 15 : (p_sigma < 40.f ? 25 : (p_sigma < 80.f ? 30 :  45))); // ASK MARC differs from manuscript
	}
	else
		o_params.nSimilarPatches = o_params.sizePatch * o_params.sizePatch * 3;

	//! Offset: step between two similar patches
	o_params.offSet = o_params.sizePatch / 2;

	//! Use the homogeneous area detection trick
	o_params.useHomogeneousArea = p_flatArea;

	//! Size of the search window around the reference patch (must be odd)
	o_params.sizeSearchWindow = o_params.nSimilarPatches / 2; // ASK MARC and 7*sizePatch1 in IPOL manuscript
	if (o_params.sizeSearchWindow % 2 == 0) o_params.sizeSearchWindow++;

	//! Search window, temporal search radii
	o_params.sizeSearchTimeRangeFwd = timeSearchRangeFwd;
	o_params.sizeSearchTimeRangeBwd = timeSearchRangeBwd;
	o_params.nSimilarPatches *= timeSearchRangeFwd + timeSearchRangeBwd + 1; // FIXME: this is just a test

	//! Size of boundaries used during the sub division
	o_params.boundary = int(1.5f * float(o_params.sizeSearchWindow));

	// TODO VIDEO: do we need to specify temporal boundary?

	//! Parameter used to determine if an area is homogeneous
	o_params.gamma = 1.05f;

	//! Parameter used to estimate the covariance matrix
	if (p_size.channels == 1)
	{
		if(s1) o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 70.f ? 1.f : 0.9f));
		else   o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 35.f ? 1.f : 0.9f));
	}
	else
	{
		if(s1) o_params.beta = 1.f;
		else   o_params.beta = (p_sigma < 50.f ? 1.2f : 1.f);
	}

	//! Parameter used to determine similar patches
	//  Differs from manuscript (tau = 4) because (1) this threshold is to be used 
	//  with the squared distance and (2) the distance is not normalized
	if(s1) o_params.tau = 0; // not used
	else   o_params.tau = 16.f * o_params.sizePatch * o_params.sizePatch * p_size.channels;

	//! Print information?
	o_params.verbose = p_verbose;

	//! Is first step?
	o_params.isFirstStep = (p_step == 1);

	//! Boost the paste trick
	o_params.doPasteBoost = true;
}

/**
 * @brief Display parameters of the NL-Bayes algorithm.
 *
 * @param i_params : nlbParams for first or second step of the algorithm;
 *
 * @return none.
 **/
void printNlbParameters(
	const nlbParams &i_prms
){
//	printf("\nVideo NLBayes parameters\n");
//	printf("------------------------\n\n");
//	printf("Noise sigma = %g\n", i_prms1.sigma);

	printf("Parameters for step %d:\n", i_prms.isFirstStep ? 1 : 2);
	printf("\tPatch search:\n");
	printf("\t\tPatch size                 = %d\n"       , i_prms.sizePatch);
	printf("\t\tNumber of patches          = %d\n"       , i_prms.nSimilarPatches);
	if (!i_prms.isFirstStep) printf("\t\tDistance threshold (tau)   = %f\n"       , i_prms.tau);
	printf("\t\tSpatial search window      = %dx%d\n"    , i_prms.sizeSearchWindow, i_prms.sizeSearchWindow);
	printf("\t\tTemporal search range      = [-%d,%d]\n" , i_prms.sizeSearchTimeRangeBwd, i_prms.sizeSearchTimeRangeBwd);
	printf("\tPatch stack filtering: beta        = %f\n" , i_prms.beta);
	printf("\tSpatial border                     = %d\n" , i_prms.boundary);
	printf("\tSpeed-ups:\n");
	printf("\t\tOffset                     = %d\n"       , i_prms.offSet);
	if (i_prms.useHomogeneousArea)
	printf("\t\tFlat area trick with gamma = %f\n"       , i_prms.gamma);
	if (i_prms.doPasteBoost)
		printf("\t\tPasteBoost (TM) active!\n\n");
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * This function splits the image in several subimages when using
 * OpenMP. Each subimage is then assinged to a thread, which runs
 * processNlBayes on the subimage. The subimages are then assembled
 * into the final image.
 *
 * @param i_noisy: contains the noisy image;
 * @param o_basic: will contain the basic estimate image after the first step;
 * @param o_final: will contain the final denoised image after the second step;
 * @param p_useArea1 : if true, use the flat area trick for the first step;
 * @param p_useArea2 : if true, use the flat area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/
int runNlBayes(
	Video_f32 const& i_noisy
,	Video_f32 &o_basic
,	Video_f32 &o_final
,	const bool p_useArea1
,	const bool p_useArea2
,	const float p_sigma
,	const bool p_verbose
){
	//! Video size
	VideoSize size = i_noisy.size();

	//! Parameters Initialization
	nlbParams p_prms1, p_prms2;
	initializeNlbParameters(p_prms1, 1, p_sigma, size, p_useArea1, p_verbose);
	initializeNlbParameters(p_prms2, 2, p_sigma, size, p_useArea2, p_verbose);

	//! NL-Bayes
	return runNlBayes(i_noisy, o_basic, o_final, p_prms1, p_prms2);
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * This function splits the image in several subimages when using
 * OpenMP. Each subimage is then assinged to a thread, which runs
 * processNlBayes on the subimage. The subimages are then assembled
 * into the final image.
 *
 * @param i_imNoisy: contains the noisy image;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_imSize: size of the image;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/
int runNlBayes(
	Video_f32 const& i_imNoisy
,	Video_f32 &o_imBasic
,	Video_f32 &o_imFinal
,	const nlbParams& p_prms1
,	const nlbParams& p_prms2
){
	//! Only 1, 3 or 4-channels images can be processed.
//	const unsigned chnls = i_imNoisy.channels;
	const unsigned chnls = 2;
	if (! (chnls == 1 || chnls == 3 || chnls == 4))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong number of channels. Must be 1, 3 or 4!!n");

	//! Number of available cores
	unsigned nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	if (p_verbose) printf("OpenMP is using %d threads\n", nThreads);
#endif
	const unsigned nParts = 2 * nThreads;

	//! Video size
	VideoSize imSize = i_imNoisy.size();

	//! Initialization
	o_imBasic.resize(imSize);
	o_imFinal.resize(imSize);

	//! Print parameters
	if (p_prms1.verbose) printNlbParameters(p_prms1);
	if (p_prms2.verbose) printNlbParameters(p_prms2);

	//! Step 1
	{
		if (p_prms1.verbose) printf("1st Step\n");

		//! RGB to YUV
		Video_f32 imNoisy = i_imNoisy;
		VideoUtils::transformColorSpace(imNoisy, true);

		//! Divide the noisy image into sub-images in order to easier parallelize the process
		// ASK MARC: any suggestion on the best way to split the space time cube
		std::vector<Video_f32> imNoisySub(nParts);
		VideoUtils::subDivide(imNoisy, imNoisySub, p_prms1.boundary, nParts);

		//! Process all sub-images
		std::vector<Video_f32> imBasicSub(nParts);
		std::vector<Video_f32> imFinalSub(nParts);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub, imSizeSub) \
		firstprivate (p_prms1)
#endif
		for (int n = 0; n < (int)nParts; n++)
			processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], p_prms1);

		//! Get the basic estimate
		VideoUtils::subBuild(imBasicSub, o_imBasic, p_prms1.boundary);

		//! YUV to RGB
		VideoUtils::transformColorSpace(o_imBasic, false);
	}

	//! Step 2
	{
		if (p_prms2.verbose) printf("2nd Step\n");

		//! Divide the noisy and basic images into sub-images in order to easier parallelize the process
		std::vector<Video_f32> imNoisySub(nParts);
		std::vector<Video_f32> imBasicSub(nParts);

		VideoUtils::subDivide(i_imNoisy, imNoisySub, p_prms2.boundary, nParts);
		VideoUtils::subDivide(o_imBasic, imBasicSub, p_prms2.boundary, nParts);

		//! Process all sub-images
		std::vector<Video_f32> imFinalSub(nParts);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub) \
		firstprivate (p_prms2)
#endif
		for (int n = 0; n < (int) nParts; n++)
			processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], p_prms2);

		//! Get the final result
		VideoUtils::subBuild(imFinalSub, o_imFinal, p_prms2.boundary);
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
void processNlBayes(
	Video_f32 const& i_imNoisy
,	Video_f32 &io_imBasic
,	Video_f32 &o_imFinal
,	nlbParams const& p_params
){
	using std::vector;

	//! Parameters initialization
	const bool step1 = p_params.isFirstStep;
	const unsigned sWx = p_params.sizeSearchWindow;
	const unsigned sWt = p_params.sizeSearchTimeRangeFwd +
	                     p_params.sizeSearchTimeRangeBwd + 1;// VIDEO
	const unsigned sP  = p_params.sizePatch;

	unsigned nInverseFailed = 0;
	const float threshold = p_params.sigma * p_params.sigma * p_params.gamma *
	                       (p_params.isFirstStep ? i_imNoisy.channels : 1.f);

	//! Weight sum per pixel // FIXME: only one channel
	Video_f32 weight(i_imNoisy.size(), 0.f);

	//! Mask: true for pixels that still need to be processed
	vector<bool> mask(i_imNoisy.wh * i_imNoisy.frames, false);

	//! Only pixels of the center of the image must be processed (not the boundaries)
	for (unsigned f =   0; f < i_imNoisy.frames      ; f++)
	for (unsigned i = sWx; i < i_imNoisy.height - sWx; i++)
	for (unsigned j = sWx; j < i_imNoisy.width  - sWx; j++)
		mask[f * i_imNoisy.wh + i * i_imNoisy.width + j] = true;

	//! Used matrices during Bayes' estimate
	const unsigned patch_dim = step1 ? sP * sP : sP * sP * i_imNoisy.channels;
	const unsigned patch_num = step1 ? p_params.nSimilarPatches : sWx * sWx * sWt;

	//! Matrices used for Bayes' estimate
	vector<unsigned> index(patch_num);
	matParams mat;
	mat.group3dTranspose.resize(patch_num * patch_dim);
	mat.tmpMat          .resize(patch_dim * patch_dim);
	mat.covMat          .resize(patch_dim * patch_dim);
	mat.covMatTmp       .resize(patch_dim * patch_dim);
	mat.baricenter      .resize(patch_dim);

	if (step1)
	{
		//! Allocate Sizes
		io_imBasic.resize(i_imNoisy.size());

		//! Matrices used for Bayes' estimate
		vector<vector<float> > group3d(i_imNoisy.channels,
				                         vector<float>(patch_num * patch_dim));

		for (unsigned ij = 0; ij < i_imNoisy.whf; ij += p_params.offSet)
			if (mask[ij]) //< Only non-seen patches are processed
			{
				if (p_params.verbose && (ij % 10000 == 0))
				{
					printf("\rprocessing step1 %05.1f", (float)ij/(float)(i_imNoisy.whf)*100.f);
					std::cout << std::flush;
				}

				const unsigned ij3 = (ij / i_imNoisy.wh) * i_imNoisy.whc + ij % i_imNoisy.wh;

				//! Search for similar patches around the reference one
				estimateSimilarPatchesStep1(i_imNoisy, group3d, index, ij3, p_params);

				//! If we use the homogeneous area trick
				bool doBayesEstimate = true;
				if (p_params.useHomogeneousArea)
					doBayesEstimate = !computeHomogeneousAreaStep1(group3d, sP, patch_num,
							threshold, i_imNoisy.size());

				//! Else, use Bayes' estimate
				if (doBayesEstimate)
					computeBayesEstimateStep1(group3d, mat, nInverseFailed, p_params);

				//! Aggregation
				computeAggregationStep1(io_imBasic, weight, mask, group3d, index, p_params);
			}

		if (p_params.verbose) printf("\n");

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, io_imBasic, weight);
	}
	else
	{
		//! Allocate Sizes
		o_imFinal.resize(i_imNoisy.size());

		//! Matrices used for Bayes' estimate
		vector<float> group3dNoisy(patch_num * patch_dim);
		vector<float> group3dBasic(patch_num * patch_dim);

		for (unsigned ij = 0; ij < i_imNoisy.whf; ij += p_params.offSet)
			if (mask[ij]) //< Only non-seen patches are processed
			{
				if (p_params.verbose && (ij % 10000 == 0))
				{
					printf("\rprocessing step2 %05.1f", (float)ij/(float)(i_imNoisy.whf)*100.f);
					std::cout << std::flush;
				}

				const unsigned ij3 = (ij / i_imNoisy.wh) * i_imNoisy.whc + ij % i_imNoisy.wh;

				//! Search for similar patches around the reference one
				unsigned nSimP = estimateSimilarPatchesStep2(i_imNoisy, io_imBasic,
						group3dNoisy, group3dBasic, index, ij3, p_params);

				//! If we use the homogeneous area trick
				bool doBayesEstimate = true;
				if (p_params.useHomogeneousArea)
					doBayesEstimate = !computeHomogeneousAreaStep2(group3dNoisy,
							group3dBasic, sP, nSimP, threshold, i_imNoisy.size());

				//! Else, use Bayes' estimate
				if (doBayesEstimate)
					computeBayesEstimateStep2(group3dNoisy, group3dBasic, mat,
							nInverseFailed, i_imNoisy.size(), p_params, nSimP);

				//! Aggregation
				computeAggregationStep2(o_imFinal, weight, mask, group3dBasic,
						index, p_params, nSimP);
			}

		if (p_params.verbose) printf("\n");

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, o_imFinal, weight);
	}

	if (nInverseFailed > 0 && p_params.verbose)
		std::cout << "nInverseFailed = " << nInverseFailed << std::endl;
}

/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param i_im: contains the noisy image on which distances are processed;
 * @param o_group3d: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void estimateSimilarPatchesStep1(
	Video_f32 const& i_im
,	std::vector<std::vector<float> > &o_group3d
,	std::vector<unsigned> &o_index
,	const unsigned pidx
,	const nlbParams &p_params
){
	//! Initialization
	const unsigned sWx   = p_params.sizeSearchWindow;
	const unsigned sWt_f = p_params.sizeSearchTimeRangeFwd;
	const unsigned sWt_b = p_params.sizeSearchTimeRangeBwd;
	const unsigned sP    = p_params.sizePatch;
	const unsigned nSimP = p_params.nSimilarPatches;

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	i_im.coords(pidx, px, py, pt, pc);

	const unsigned rangex[2] = {px - (sWx - 1)/2, px + (sWx - 1)/2};
	const unsigned rangey[2] = {py - (sWx - 1)/2, py + (sWx - 1)/2};
	const unsigned ranget[2] = {std::max((int)pt - (int)sWt_b, 0),
	                            std::min((int)pt + (int)sWt_f, (int)i_im.frames-1)};

	std::vector<std::pair<float, unsigned> > distance(sWx * sWx * 
	                                            (ranget[1] - ranget[0] + 1));

	//! Compute distance between patches in search range
	for (unsigned qt = ranget[0], dt = 0; qt <= ranget[1]; qt++, dt++)
	for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
	for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
	{
		//! Squared L2 distance
		float diff = 0.f;
		for (unsigned hy = 0; hy < sP; hy++)
		for (unsigned hx = 0; hx < sP; hx++)
		{
			const float tmp = i_im(px + hx, py + hy, pt)
			                - i_im(qx + hx, qy + hy, qt);
			diff += tmp * tmp;
		}

		//! Save distance and corresponding patch index
		distance[dt * sWx*sWx + dy * sWx + dx] = 
			std::make_pair(diff, i_im.index(qx, qy, qt, 0));
	}

	//! Keep only the N2 best similar patches
	std::partial_sort(distance.begin(), distance.begin() + nSimP,
	                  distance.end(), comparaisonFirst);

	//! Register position of patches
	for (unsigned n = 0; n < nSimP; n++) o_index[n] = distance[n].second;

	//! Stack selected patches into the 3D group
	const unsigned w  = i_im.width;
	const unsigned wh = i_im.wh;
	for (unsigned c  = 0; c < i_im.channels; c++)
	for (unsigned hy = 0, k = 0; hy < sP; hy++)
	for (unsigned hx = 0;        hx < sP; hx++)
	for (unsigned n  = 0; n < nSimP; n++, k++)
		o_group3d[c][k] = i_im.data[o_index[n] + hy * w + hx + c * wh];

	/* 00  pixels from all patches
	 * 01  pixels from all patches
	 * ...
	 * 0sp pixels from all patches
	 */
}

/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param i_imBasic: contains the basic estimation;
 * @param o_group3dNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_group3dBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_imSize: size of images;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video_f32 const& i_imNoisy
,	Video_f32 const& i_imBasic
,	std::vector<float> &o_group3dNoisy
,	std::vector<float> &o_group3dBasic
,	std::vector<unsigned> &o_index
,	const unsigned pidx
,	const nlbParams &p_params
){
	//! Initialization
	const unsigned width = i_imNoisy.width;
	const unsigned chnls = i_imNoisy.channels;
	const unsigned wh    = width * i_imNoisy.height;
	const unsigned sP    = p_params.sizePatch;
	const unsigned sWx   = p_params.sizeSearchWindow;
	const unsigned sWt_f = p_params.sizeSearchTimeRangeFwd;
	const unsigned sWt_b = p_params.sizeSearchTimeRangeBwd;

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	i_imBasic.coords(pidx, px, py, pt, pc);

	const unsigned rangex[2] = {px - (sWx - 1)/2, px + (sWx - 1)/2};
	const unsigned rangey[2] = {py - (sWx - 1)/2, py + (sWx - 1)/2};
	const unsigned ranget[2] = {std::max((int)pt - (int)sWt_b, 0),
	                            std::min((int)pt + (int)sWt_f, (int)i_imNoisy.frames-1)};

	std::vector<std::pair<float, unsigned> > distance(sWx * sWx * 
	                                            (ranget[1] - ranget[0] + 1));

	//! Compute distance between patches in search range
	for (unsigned qt = ranget[0], dt = 0; qt <= ranget[1]; qt++, dt++)
	for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
	for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
	{
		//! Squared L2 distance
		float diff = 0.f;
		for (unsigned c = 0; c < chnls; c++)
		for (unsigned hy = 0; hy < sP; hy++)
		for (unsigned hx = 0; hx < sP; hx++)
		{
			const float tmp = i_imBasic(px + hx, py + hy, pt, c)
			                - i_imBasic(qx + hx, qy + hy, qt, c);
			diff += tmp * tmp;
		}

		//! Save distance and corresponding patch index
		distance[dt * sWx*sWx + dy * sWx + dx] = 
			std::make_pair(diff, i_imBasic.index(qx, qy, qt, 0));
	}

	//! Keep only the nSimilarPatches best similar patches
	std::partial_sort(distance.begin(), distance.begin() + p_params.nSimilarPatches,
	                  distance.end(), comparaisonFirst);

	//! Save index of similar patches
	const float threshold = (p_params.tau > distance[p_params.nSimilarPatches - 1].first ?
	                         p_params.tau : distance[p_params.nSimilarPatches - 1].first);
	unsigned nSimP = 0;

	//! Register position of similar patches
	for (unsigned n = 0; n < distance.size(); n++)
		if (distance[n].first < threshold)
			o_index[nSimP++] = distance[n].second;

	//! Save similar patches into 3D groups
	for (unsigned c = 0, k = 0; c < chnls; c++)
	for (unsigned hy = 0; hy < sP; hy++)
	for (unsigned hx = 0; hx < sP; hx++)
	for (unsigned n = 0; n < nSimP; n++, k++)
	{
		o_group3dNoisy[k] = i_imNoisy.data[c * wh + o_index[n] + hy * width + hx];
		o_group3dBasic[k] = i_imBasic.data[c * wh + o_index[n] + hy * width + hx];
	}

	return nSimP;
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3d: contains for each channels values of similar patches. If an homogeneous area
 *			is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_doLinearRegression: if true, apply a linear regression to average value of pixels;
 * @param p_imSize: size of the image.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep1(
	std::vector<std::vector<float> > &io_group3d
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
){
	//! Initialization
	const unsigned N = p_sP * p_sP * p_nSimP;

	//! Compute the standard deviation of the set of patches
	float stdDev = 0.f;
	for (unsigned c = 0; c < p_imSize.channels; c++)
		stdDev += computeStdDeviation(io_group3d[c], p_sP * p_sP, p_nSimP, 1);

	//! If we are in an homogeneous area
	if (stdDev < p_threshold)
	{
		for (unsigned c = 0; c < p_imSize.channels; c++)
		{
			float mean = 0.f;

			for (unsigned k = 0; k < N; k++) mean += io_group3d[c][k];
			mean /= (float) N;

			for (unsigned k = 0; k < N; k++) io_group3d[c][k] = mean;
		}
		return 1;
	}
	else return 0;
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3dNoisy: contains values of similar patches for the noisy image;
 * @param io_group3dBasic: contains values of similar patches for the basic image. If an homogeneous
 *		area is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the image.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep2(
	std::vector<float> const& i_group3dNoisy
,	std::vector<float> &io_group3dBasic
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
){
	//! Parameters
	const unsigned sP2 = p_sP * p_sP;
	const unsigned sPC = sP2 * p_imSize.channels;

	//! Compute the standard deviation of the set of patches
	const float stdDev = computeStdDeviation(i_group3dNoisy, sP2, p_nSimP, p_imSize.channels);

	//! If we are in an homogeneous area
	if (stdDev < p_threshold)
	{
		for (unsigned c = 0; c < p_imSize.channels; c++)
		{
				float mean = 0.f;
				for (unsigned n = 0; n < p_nSimP; n++)
				for (unsigned k = 0; k < sP2; k++)
					mean += io_group3dBasic[n * sPC + c * sP2 + k];

				mean /= float(sP2 * p_nSimP);

				for (unsigned n = 0; n < p_nSimP; n++)
				for (unsigned k = 0; k < sP2; k++)
					io_group3dBasic[n * sPC + c * sP2 + k] = mean;
		}
		return 1;
	}
	else return 0;
}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_group3d: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *    - group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void computeBayesEstimateStep1(
	std::vector<std::vector<float> > &io_group3d
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
){
	//! Parameters
	const unsigned chnls = io_group3d.size();
	const unsigned nSimP = p_params.nSimilarPatches;
	const unsigned sP2   = p_params.sizePatch * p_params.sizePatch;
	const float valDiag  = p_params.beta * p_params.sigma * p_params.sigma;

	//! Bayes estimate
	for (unsigned c = 0; c < chnls; c++)
	{
		//! Center data around the baricenter
		centerData(io_group3d[c], i_mat.baricenter, nSimP, sP2);

		//! Compute the covariance matrix of the set of similar patches
		covarianceMatrix(io_group3d[c], i_mat.covMat, nSimP, sP2);

		//! Bayes' Filtering
		if (inverseMatrix(i_mat.covMat, sP2) == EXIT_SUCCESS)
		{
			productMatrix(i_mat.group3dTranspose, i_mat.covMat, io_group3d[c],
			              sP2, sP2, nSimP);
			for (unsigned k = 0; k < sP2 * nSimP; k++)
				io_group3d[c][k] -= valDiag * i_mat.group3dTranspose[k];
		}
		else 
			io_nInverseFailed++;

		//! Add baricenter
		for (unsigned j = 0, k = 0; j < sP2; j++)
			for (unsigned i = 0; i < nSimP; i++, k++)
				io_group3d[c][k] += i_mat.baricenter[j];
	}
}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param i_group3dNoisy: contains all similar patches in the noisy image;
 * @param io_group3dBasic: contains all similar patches in the basic image. Will contain estimates
 *       for all similar patches;
 * @param i_mat: contains :
 *    - group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_imSize: size of the image;
 * @param p_params: see processStep2 for more explanations;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeBayesEstimateStep2(
	std::vector<float> &i_group3dNoisy
,	std::vector<float> &io_group3dBasic
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float diagVal = p_params.beta * p_params.sigma * p_params.sigma;
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch * p_imSize.channels;

	//! Center 3D groups around their baricenter
	centerData(io_group3dBasic, i_mat.baricenter, p_nSimP, sPC);
	centerData(i_group3dNoisy, i_mat.baricenter, p_nSimP, sPC);

	//! Compute the covariance matrix of the set of similar patches
	covarianceMatrix(io_group3dBasic, i_mat.covMat, p_nSimP, sPC);

	//! Bayes' Filtering
	for (unsigned k = 0; k < sPC; k++)
		i_mat.covMat[k * sPC + k] += diagVal;

	//! Compute the estimate
	if (inverseMatrix(i_mat.covMat, sPC) == EXIT_SUCCESS)
	{
		productMatrix(io_group3dBasic, i_mat.covMat, i_group3dNoisy, sPC, sPC, p_nSimP);
		for (unsigned k = 0; k < sPC * p_nSimP; k++)
			io_group3dBasic[k] = i_group3dNoisy[k] - diagVal * io_group3dBasic[k];
	}
	else 
		io_nInverseFailed++;

	//! Add baricenter
	for (unsigned j = 0, k = 0; j < sPC; j++)
		for (unsigned i = 0; i < p_nSimP; i++, k++)
			io_group3dBasic[k] += i_mat.baricenter[j];
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void computeAggregationStep1(
	Video_f32 &io_im
,	Video_f32 &io_weight
,	std::vector<bool> &io_mask
,	std::vector<std::vector<float> > const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
){
	//! Parameters initializations
	const unsigned chnls  = io_im.channels;
	const unsigned width  = io_im.width;
	const unsigned height = io_im.height;
	const unsigned sP     = p_params.sizePatch;
	const unsigned nSimP  = p_params.nSimilarPatches;

	//! Aggregate estimates
	for (unsigned n = 0; n < nSimP; n++)
	{
		const unsigned ind = i_index[n];
		for (unsigned c = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * width * height;
			for (unsigned p = 0; p < sP; p++)
			for (unsigned q = 0; q < sP; q++)
			{
				io_im.data[ij + p * width + q] += i_group3d[c][(p * sP + q) * nSimP + n];
				io_weight.data[ij + p * width + q]++;
			}
		}

		//! Use Paste Trick
		const unsigned ind1 = (ind / io_im.whc) * io_im.wh + ind % io_im.wh;
		io_mask[ind1] = false;

		if (p_params.doPasteBoost)
		{
			io_mask[ind1 - width] = false;
			io_mask[ind1 + width] = false;
			io_mask[ind1 - 1    ] = false;
			io_mask[ind1 + 1    ] = false;
		}
	}
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeAggregationStep2(
	Video_f32 &io_im
,	Video_f32 &io_weight
,	std::vector<bool> &io_mask
,	std::vector<float> const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls = io_im.channels;
	const unsigned width = io_im.width;
	const unsigned wh    = width * io_im.height;
	const unsigned sP    = p_params.sizePatch;

	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned p = 0; p < sP; p++)
			for (unsigned q = 0; q < sP; q++, k++)
			{
				io_im.data[ij + p * width + q] += i_group3d[k * p_nSimP + n];
				io_weight.data[ij + p * width + q]++;
			}
		}

		//! Apply Paste Trick
		const unsigned ind1 = (ind / io_im.whc) * io_im.wh + ind % io_im.wh;
		io_mask[ind1] = false;

		if (p_params.doPasteBoost)
		{
			io_mask[ind1 - width ] = false;
			io_mask[ind1 + width ] = false;
			io_mask[ind1 - 1     ] = false;
			io_mask[ind1 + 1     ] = false;
		}
	}
}

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_im: image of reference, when the weight if null;
 * io_im: will contain the final image;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return : none.
 **/
void computeWeightedAggregation(
	Video_f32 const& i_im
,	Video_f32 &io_im
,	Video_f32 const& i_weight
){
	for (unsigned k = 0; k < i_im.whcf; k++)
		//! To avoid weighting problem (particularly near boundaries of the image)
		if (i_weight(k) > 0.f) io_im(k) /= i_weight(k);
		else                   io_im(k) /= i_im(k);

}

}
