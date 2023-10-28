# Active_Finite_Voronoi_simulation

This is the core simulation MATLAB code of the [active finite Voronoi paper]
This repo includes the core MATLAB codes for our paper: ["Bridging the Gap Between Collective Motility and Epithelial-Mesenchymal Transitions through the Active Finite Voronoi Model", Junxiang Huang, Herbert Levine, Dapeng Bi](https://pubs.rsc.org/en/content/articlelanding/2023/sm/d3sm00327b).

Some sample movies can be viewed here: https://drive.google.com/drive/folders/1hlINXSVOMjMYvTwKG0yyWhYWvtaviSAn?usp=share_link

The source code files included in this repo are:
* AFV_code_illustration.mlx: The main demo notebook/live script. One can easily run our code interactively and check visualization results. 
* make_finite_voronoi_pbc.m: supportive function to generate finite Voronoi tissue
* get_finite_voronoi_force.m: calculate the interaction force (not include v0)
* draw_finite_voronoi.m: visualize the cells in a PBC box
  
To run the code, please download everything and put them in the same folder, and add the folder to your MATLAB path. Then open “AFV_code_illustration.mlx” and click the “run to end” button in your MATLAB interface. The resulting output is expected to be the same (a pdf archived version of the notebook is also provided for easy comparison).

For any questions regarding the code, please contact Junxiang Huang (jxhuangphys@gmail.com).

You are also welcome to contact the PI's for shared interests:
* Max Bi (email: d.bi@northeastern.edu, homepage: https://sites.google.com/view/dapengbi)
* Herbert Levine (email: h.levine@northeastern.edu, profile: https://scholar.google.com/citations?user=gHjHIvAAAAAJ).
