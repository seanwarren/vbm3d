/*
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work: Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <mex.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <string>
#include <sstream>

#include "utils/utilities.h"
#include "nldct.h"
#include "utils/cmd_option.h"

using namespace std;

enum Mode { BSIC_DENO, BSIC_ONLY, DENO_ONLY, NISY_ONLY };

class clo
{
public:
   clo(int nrhs, const mxArray *prhs[]) : 
      nrhs(nrhs), prhs(prhs)
   {
   }

   const mxArray* get(const char *o)
   {
      const int buflen = 1024;
      char v[buflen];

      for (int i = 0; i < nrhs - 1; i++)
         if (mxIsChar(prhs[i]))
         {
            mxGetString(prhs[i], v, buflen);
            if (0 == strcmp(v, o)) return prhs[i + 1];
         }

      return nullptr;
   }

   bool has(const char *o)
   {
      const int buflen = 1024;
      char v[buflen];

      for (int i = 0; i < nrhs; i++)
         if (mxIsChar(prhs[i]))
         {
            mxGetString(prhs[i], v, buflen);
            if (0 == strcmp(v, o)) return true;
         }

      return false;
   }

   template<typename T>
   T pick(const char* o, T default, const char* desc = nullptr)
   {
      const mxArray* v = get(o);
      if (v) return (T)mxGetScalar(v);
      return default;
   }

private:
   int nrhs;
   const mxArray **prhs;
};


Video<float> videoFromMatlab(const mxArray* v)
{
   mxAssert(mxIsSingle(v), "Input image must be of type single");
   int ndim = mxGetNumberOfDimensions(v);
   const mwSize* dims = mxGetDimensions(v);
   int width = dims[0];
   int height = dims[1];

   int chnls = 1;
   int frames = 1;

   if (ndim >= 3)
      frames = dims[2];
   if (ndim >= 4)
      chnls = dims[3];

   Video<float> vid(width, height, frames, chnls);
   float* data = (float*)mxGetData(v);
   std::copy_n(data, width*height*frames*chnls, vid.data.begin());
   return vid;
}

mxArray* videoToMatlab(const Video<float>& vid)
{
   auto s = vid.sz;
   mwSize dims[4] = { s.width, s.height, s.frames, s.channels };
   mxArray* v = mxCreateNumericArray(4, dims, mxSINGLE_CLASS, mxREAL);
   std::copy(vid.data.begin(), vid.data.end(), (float*)mxGetData(v));
   return v;
}


void mexFunction(int nlhs, mxArray *plhs[],
   int nrhs, const mxArray *prhs[])
{
   /*
	clo_usage("Video NL-Bayes video denoising");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");
   */

   clo c(nrhs, prhs);

   /*
	//! Paths to input/output sequences
	using std::string;
	const string  input_path = c.pick("i"    , ""              , "< input sequence");
	const string  inbsc_path = c.pick("b"    , ""              , "< input basic sequence");
	const string  noisy_path = c.pick("nisy" , ""              , "> noisy sequence");
	const string  final_path = c.pick("deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = c.pick("bsic" , "bsic_%03d.png" , "> basic denoised sequence");
	const string   diff_path = c.pick("diff" , "diff_%03d.png" , "> difference sequence");

	const unsigned firstFrame = c.pick("f", 0, "first frame");
	const unsigned lastFrame  = c.pick("l", 0, "last frame");
	const unsigned frameStep  = c.pick("s", 1, "frame step");
   */

	//! General parameters
	const float sigma = c.pick("sigma", 0.f, "Add noise of standard deviation sigma");
   const bool has_noise = c.has("has-noise");
   const bool verbose = c.has("verbose");
   const bool order_inv = c.has("order-inv"); // "> invariance to patch order in stack"
	const unsigned print_prms = (unsigned) c.pick("print-prms", 0, "> prints parameters for given channels");

	//! Video NLB parameters
	const int time_search1  = c.pick("wt1", 0  , "> Search window temporal radius, step 1");
	const int time_search2  = c.pick("wt2", 0  , "> Search window temporal radius, step 2");
	const int space_search1 = c.pick("wx1",-1  , "> Search window spatial radius, step 1");
	const int space_search2 = c.pick("wx2",-1  , "> Search window spatial radius, step 2");
	const int patch_sizex1  = c.pick("px1",-1  , "> Spatial patch size, step 1");
	const int patch_sizex2  = c.pick("px2",-1  , "> Spatial patch size, step 2");
	const int patch_sizet1  = c.pick("pt1", 1  , "> Temporal patch size, step 1");
	const int patch_sizet2  = c.pick("pt2", 1  , "> Temporal patch size, step 2");
	const int num_patches1  = c.pick("np1",-1  , "> Number of similar patches, step 1");
	const int num_patches2  = c.pick("np2",-1  , "> Number of similar patches, step 2");
	const float tau1        = c.pick("t1" ,-1.f, "> Step 1 distance threshold");
	const float tau2        = c.pick("t2" ,-1.f, "> Step 2 distance threshold");
	const float beta1       = c.pick("b1" ,-1.f, "> Noise correction factor beta, step 1");
	const float beta2       = c.pick("b2" ,-1.f, "> Noise correction factor beta, step 2");
#ifdef VBM3D_SEARCH
	const int space_search_f1= c.pick("wxf1",-1  , "> Search window for predictive search, step 1");
	const int space_search_f2= c.pick("wxf2",-1  , "> Search window for predictive search, step 2");
	const int num_patches_f1 = c.pick("npf1",-1  , "> Number of similar patches per frame, step 1");
	const int num_patches_f2 = c.pick("npf2",-1  , "> Number of similar patches per frame, step 2");
	const float dsub1        = c.pick("dsub1" ,-1.f, "> Distance bias, step 1");
	const float dsub2        = c.pick("dsub2" ,-1.f, "> Distance bias, step 2");
#endif
	const float beta_mean1  = c.pick("bm1",-1.f, "> Noise correction factor beta for mean, step 1");
	const float beta_mean2  = c.pick("bm2",-1.f, "> Noise correction factor beta for mean, step 2");
	const bool no_paste1  = c.pick("no-paste1", false , "> disable paste trick, step 1");
	const bool no_paste2  = c.pick("no-paste2", false , "> disable paste trick, step 2");
	const bool no_step1   = c.pick("no-step1" , false , "> disable patch skipping, step 1");
	const bool no_step2   = c.pick("no-step2" , false , "> disable patch skipping, step 2");
	const bool agg_win1   = c.pick("agg-win1" , false , "> aggregation window, step 1");
	const bool agg_win2   = c.pick("agg-win2" , false , "> aggregation window, step 2");


	if ((patch_sizex1 <= 0 && patch_sizex1 != -1) ||
	    (patch_sizex2 <= 0 && patch_sizex2 != -1) ||
	    (patch_sizet1 <= 0 || patch_sizet2 <=   0) )
	{
      mexWarnMsgIdAndTxt("vnldct:error", "px1, px2, pt1 and pt2 cannot be negative.");
		return;
	}

	if ((num_patches1 < 0 && num_patches1 != -1) ||
	    (num_patches2 < 0 && num_patches2 != -1) )
	{
      mexWarnMsgIdAndTxt("vnldct:error", "np1 and np2 cannot be negative.");
		return;
	}

	if ((space_search1 < 0 && space_search1 != -1) ||
	    (space_search2 < 0 && space_search2 != -1) ||
	    ( time_search1 < 0 ||  time_search2 <   0) )
	{
      mexWarnMsgIdAndTxt("vnldct:error", "wx1, wx2, wt1 and wt2 cannot be negative.");
		return;
	}

   const mxArray* bflow_mat = c.get("bflow");
   const mxArray* fflow_mat = c.get("fflow");

	if ((fflow_mat && !bflow_mat) || (!fflow_mat && bflow_mat))
	{
      mexWarnMsgIdAndTxt("vnldct:error", "Only one oflow path provided.");
		return;
	}

	//! Only print parameters
	if (print_prms)
	{
		VideoSize tmp;
		tmp.channels = print_prms;

		//! Compute denoising default parameters
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, sigma, tmp, verbose);
		VideoNLB::initializeNlbParameters(prms2, 2, sigma, tmp, verbose);

		//! Override with command line parameters
		if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
		if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
		if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, tmp, (unsigned)patch_sizex1);;
		if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, tmp, (unsigned)patch_sizex2);;
		if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
		if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
		if (beta1         >= 0) prms1.beta = beta1;
		if (beta2         >= 0) prms2.beta = beta2;
		if (beta1         >= 0) prms1.betaMean = beta1;
		if (beta2         >= 0) prms2.betaMean = beta2;
		if (beta_mean1    >= 0) prms1.betaMean = beta_mean1;
		if (beta_mean2    >= 0) prms2.betaMean = beta_mean2;

		if (patch_sizet1  >= 0) prms1.sizePatchTime = patch_sizet1;
		if (patch_sizet2  >= 0) prms2.sizePatchTime = patch_sizet2;

#ifdef VBM3D_SEARCH
		if (space_search_f1 >= 0) prms1.sizeSearchWindowPred = space_search_f1;
		if (space_search_f2 >= 0) prms2.sizeSearchWindowPred = space_search_f2;
		if (num_patches_f1  >= 0) prms1.nSimilarPatchesPred = num_patches_f1;
		if (num_patches_f2  >= 0) prms2.nSimilarPatchesPred = num_patches_f2;
		if (dsub1           >= 0) prms1.dsub = dsub1;
		if (dsub2           >= 0) prms2.dsub = dsub2;
		if (tau1            >= 0) prms1.tau = tau1;
		if (tau2            >= 0) prms2.tau = tau2;
#endif
		if (agg_win1) prms1.agg_window = true;
		if (agg_win2) prms2.agg_window = true;

		if (no_paste1) prms1.doPasteBoost = false;
		if (no_paste2) prms2.doPasteBoost = false;
		if (no_step1)  prms1.offSet = 1;
		if (no_step2)  prms2.offSet = 1;

		if (order_inv) prms1.orderInvariance = true;
		if (order_inv) prms2.orderInvariance = true;
 
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		return;
	}

	//! Declarations
	Video<float> noisy, basic, final;
	Video<float> fflow, bflow;

	//! Load input videos
   Video<float> original = videoFromMatlab(prhs[0]);
   if (fflow_mat) fflow = videoFromMatlab(fflow_mat);
   if (bflow_mat) bflow = videoFromMatlab(bflow_mat);

	//! Denoising
	if (verbose) printf("Running Video NL-Bayes on the noisy video\n");

	//! Compute denoising default parameters
	VideoNLB::nlbParams prms1, prms2;
	VideoNLB::initializeNlbParameters(prms1, 1, sigma, noisy.sz, verbose);
	VideoNLB::initializeNlbParameters(prms2, 2, sigma, noisy.sz, verbose);

	//! Override with command line parameters
	if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
	if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
	if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, noisy.sz, (unsigned)patch_sizex1);;
	if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, noisy.sz, (unsigned)patch_sizex2);;
	if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
	if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
	if (beta1         >= 0) prms1.beta = beta1;
	if (beta2         >= 0) prms2.beta = beta2;
	if (beta1         >= 0) prms1.betaMean = beta1;
	if (beta2         >= 0) prms2.betaMean = beta2;
	if (beta_mean1    >= 0) prms1.betaMean = beta_mean1;
	if (beta_mean2    >= 0) prms2.betaMean = beta_mean2;

	if (patch_sizet1  >= 0) prms1.sizePatchTime = patch_sizet1;
	if (patch_sizet2  >= 0) prms2.sizePatchTime = patch_sizet2;

#ifdef VBM3D_SEARCH
	if (space_search_f1 >= 0) prms1.sizeSearchWindowPred = space_search_f1;
	if (space_search_f2 >= 0) prms2.sizeSearchWindowPred = space_search_f2;
	if (num_patches_f1  >= 0) prms1.nSimilarPatchesPred = num_patches_f1;
	if (num_patches_f2  >= 0) prms2.nSimilarPatchesPred = num_patches_f2;
	if (dsub1           >= 0) prms1.dsub = dsub1;
	if (dsub2           >= 0) prms2.dsub = dsub2;
	if (tau1            >= 0) prms1.tau = tau1;
	if (tau2            >= 0) prms2.tau = tau2;
#endif
	if (agg_win1) prms1.agg_window = true;
	if (agg_win2) prms2.agg_window = true;

	if (no_paste1) prms1.doPasteBoost = false;
	if (no_paste2) prms2.doPasteBoost = false;
	if (no_step1)  prms1.offSet = 1;
	if (no_step2)  prms2.offSet = 1;

	if (order_inv) prms1.orderInvariance = true;
	if (order_inv) prms2.orderInvariance = true;

	//! Percentage or processed groups of patches over total number of pixels
	std::vector<float> groupsRatio;

	//! Run denoising algorithm
#ifndef DEBUG_COMPUTE_GROUP_ERROR
	if (use_oflow)
		groupsRatio = VideoNLB::runNlBayes(noisy, fflow, bflow, basic, final, prms1, prms2);
	else
		groupsRatio = VideoNLB::runNlBayes(noisy, basic, final, prms1, prms2);
#else
	if (use_oflow)
		groupsRatio = VideoNLB::runNlBayes(noisy, fflow, bflow, basic, final, prms1, prms2, original);
	else
		groupsRatio = VideoNLB::runNlBayes(noisy, basic, final, prms1, prms2, original);
#endif

   if (nlhs >= 0)
      plhs[0] = videoToMatlab(final);

	if (verbose) printf("Done\n");
	return;
}
