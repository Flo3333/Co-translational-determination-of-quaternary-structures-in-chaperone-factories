import CustomPandasFramework.FishBicolor.analysis as analysis


input_path = '/media/floricslimani/SSD 4To/SSD_floricslimani/Fish_Bicouleur/analysis_input/'
output_path = input_path.replace('input', 'output')

do_colocalisation_extract = True

do_clustered_spots_extract = True

do_colocalisation_plots = True
colocalisation_measure_list = ['rna1_number', 'rna2_number', 'suntag_number', 'rna1_colocalising_with_rna2_fraction', 'rna1_colocalising_with_suntag_fraction', 'rna2_colocalising_with_rna1_fraction', 'rna2_colocalising_with_rna1_fraction', 'rna2_colocalising_with_suntag_fraction', 'rna1_closer_1000nm_rna2_fraction', 'rna1_closer_2000nm_rna2_fraction']

do_clustered_spots_plots = True

do_distance_closest_suntag = True

#Analysis
print("Starting analysis")
if do_colocalisation_extract :
    print('\ncolocalisation excel extract...')
    analysis.colocalisation_extract(input_path, output_path)

if do_clustered_spots_extract :
    print('\nclustered_spots excel extract...')
    analysis.cluster_quantification_extract(input_path, output_path)

if do_colocalisation_plots :
    print('\nStarting colocalisation_plots :')
    analysis.colocalisation_plots(input_path, output_path, measure_list= colocalisation_measure_list)

if do_clustered_spots_plots :
    print('\nStarting clusters_plot :')
    analysis.clusters_plots(input_path, output_path, measure_list= ['free_spot_number','mean_suntag_per_free_spot','std_suntag_per_free_spot','median_suntag_per_free_spot','total_suntag_colocalising_free_spot','clustered_spot_number','mean_suntag_per_clusterd_spot','std_suntag_per_clusterd_spot','median_suntag_per_clusterd_spot','total_suntag_colocalising_clusterd_spot'])

if do_distance_closest_suntag :
    print("Distance to closest suntag distribution...")
    analysis.closest_suntag_distance(input_path, output_path)