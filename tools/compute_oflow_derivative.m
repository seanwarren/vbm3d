% computes convective derivative of a result on a sequence of sintel database

% add folders to path
addpath ~/Work/denoising/projects/video_nlbayes3d/tools/
addpath ~/Work/doctorado/codigo/algos/video_poisson/rida/vp_clean/matlab/paper_code
addpath ~/Work/doctorado/codigo/algos/optical_flow/flow-code-matlab/               

% sequence name
seq = 'bandage_1';
sigma = '40';
pt = '3';

root_folder = '/media/pariasm/tera/funes/denoising/';

% flow data folder
orig_pat = [root_folder 'data/sintel_training/final/'      seq '/frame_%04d.png'];
occl_pat = [root_folder 'data/sintel_training/occlusions/' seq '/frame_%04d.png'];
flow_pat = [root_folder 'data/sintel_training/flow/'       seq '/frame_%04d.flo'];

% bm4d pattern 
bm4d_pat = [root_folder 'algos/results/VBM4D/' seq '_mp_s' sigma '/deno_%03d.png'];
vnlb_pat = [root_folder 'projects/video_nlbayes3d/results/vnlbayes/' ... 
                        'table_1ps5_2ps5_2wx37_2r16_2np160_1r16_1np15sigma/' ...
                        seq '_s' sigma '_pt' pt '/deno_%03d.png'];

% load denoised sequence
orig =    load_sequence(orig_pat,1,50);
deno =    load_sequence(vnlb_pat,1,50);
occl =    load_sequence(occl_pat,1,49);
flow = readFlowSequence(flow_pat,1,50);

% diff = deno - orig;
%diff = orig;
diff = deno;

% compute convective derivative
flow_x = squeeze(flow(:,:,1,:));
flow_y = squeeze(flow(:,:,2,:));
occl = 1 - occl/255;

%Dflow = convective_derivative(flow_x, flow_y, [], [], occl, [], 'forward');   

% filter sequences
for ch = 1:size(diff,3),
	chnl = double(diff(:,:,ch,1:49));    
	chnl = single(Dflow * chnl(:));
	diff(:,:,ch,1:49) = reshape(chnl, size(diff(:,:,ch,1:49)));
end

round(20*log10(255/sqrt(mean(diff(:).^2)))*100)/100
