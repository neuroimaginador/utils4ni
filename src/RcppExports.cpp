// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// adaptive_nonlocal_means_denoising
NumericVector adaptive_nonlocal_means_denoising(NumericVector ima, IntegerVector dims, int v, int f, double h);
RcppExport SEXP _utils4ni_adaptive_nonlocal_means_denoising(SEXP imaSEXP, SEXP dimsSEXP, SEXP vSEXP, SEXP fSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ima(imaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type f(fSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(adaptive_nonlocal_means_denoising(ima, dims, v, f, h));
    return rcpp_result_gen;
END_RCPP
}
// cCreateMask
NumericVector cCreateMask(NumericVector ima, IntegerVector dims, int c, NumericVector medias);
RcppExport SEXP _utils4ni_cCreateMask(SEXP imaSEXP, SEXP dimsSEXP, SEXP cSEXP, SEXP mediasSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ima(imaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type c(cSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type medias(mediasSEXP);
    rcpp_result_gen = Rcpp::wrap(cCreateMask(ima, dims, c, medias));
    return rcpp_result_gen;
END_RCPP
}
// cgradientmodule
NumericVector cgradientmodule(NumericVector ima, IntegerVector dims);
RcppExport SEXP _utils4ni_cgradientmodule(SEXP imaSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ima(imaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(cgradientmodule(ima, dims));
    return rcpp_result_gen;
END_RCPP
}
// connected_components
IntegerVector connected_components(IntegerVector image);
RcppExport SEXP _utils4ni_connected_components(SEXP imageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type image(imageSEXP);
    rcpp_result_gen = Rcpp::wrap(connected_components(image));
    return rcpp_result_gen;
END_RCPP
}
// cTruncado
NumericVector cTruncado(NumericVector ima, IntegerVector dims);
RcppExport SEXP _utils4ni_cTruncado(SEXP imaSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ima(imaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(cTruncado(ima, dims));
    return rcpp_result_gen;
END_RCPP
}
// cvolfilter2d
NumericVector cvolfilter2d(NumericVector ima, IntegerVector dims, int v);
RcppExport SEXP _utils4ni_cvolfilter2d(SEXP imaSEXP, SEXP dimsSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ima(imaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(cvolfilter2d(ima, dims, v));
    return rcpp_result_gen;
END_RCPP
}
// defuzzify
IntegerVector defuzzify(NumericVector image);
RcppExport SEXP _utils4ni_defuzzify(SEXP imageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type image(imageSEXP);
    rcpp_result_gen = Rcpp::wrap(defuzzify(image));
    return rcpp_result_gen;
END_RCPP
}
// extend_labels
IntegerVector extend_labels(IntegerVector pIn, IntegerVector maskImage);
RcppExport SEXP _utils4ni_extend_labels(SEXP pInSEXP, SEXP maskImageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type pIn(pInSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type maskImage(maskImageSEXP);
    rcpp_result_gen = Rcpp::wrap(extend_labels(pIn, maskImage));
    return rcpp_result_gen;
END_RCPP
}
// map_ids_workhorse
IntegerVector map_ids_workhorse(IntegerVector x, IntegerVector source, IntegerVector target);
RcppExport SEXP _utils4ni_map_ids_workhorse(SEXP xSEXP, SEXP sourceSEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type source(sourceSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(map_ids_workhorse(x, source, target));
    return rcpp_result_gen;
END_RCPP
}
// map_extra_classes
IntegerVector map_extra_classes(IntegerVector x, IntegerVector source, int remaining);
RcppExport SEXP _utils4ni_map_extra_classes(SEXP xSEXP, SEXP sourceSEXP, SEXP remainingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type source(sourceSEXP);
    Rcpp::traits::input_parameter< int >::type remaining(remainingSEXP);
    rcpp_result_gen = Rcpp::wrap(map_extra_classes(x, source, remaining));
    return rcpp_result_gen;
END_RCPP
}
// mask_values
NumericVector mask_values(NumericVector input, double low, double high, double in_val, double out_val);
RcppExport SEXP _utils4ni_mask_values(SEXP inputSEXP, SEXP lowSEXP, SEXP highSEXP, SEXP in_valSEXP, SEXP out_valSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input(inputSEXP);
    Rcpp::traits::input_parameter< double >::type low(lowSEXP);
    Rcpp::traits::input_parameter< double >::type high(highSEXP);
    Rcpp::traits::input_parameter< double >::type in_val(in_valSEXP);
    Rcpp::traits::input_parameter< double >::type out_val(out_valSEXP);
    rcpp_result_gen = Rcpp::wrap(mask_values(input, low, high, in_val, out_val));
    return rcpp_result_gen;
END_RCPP
}
// get_neighbours
IntegerVector get_neighbours(NumericVector array, int width);
RcppExport SEXP _utils4ni_get_neighbours(SEXP arraySEXP, SEXP widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type array(arraySEXP);
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    rcpp_result_gen = Rcpp::wrap(get_neighbours(array, width));
    return rcpp_result_gen;
END_RCPP
}
// fast_read_nifti
NumericVector fast_read_nifti(SEXP filename);
RcppExport SEXP _utils4ni_fast_read_nifti(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_read_nifti(filename));
    return rcpp_result_gen;
END_RCPP
}
// count_elegible
int count_elegible(NumericVector image, int patch_size, int search_size, int stride, IntegerVector voxel_lookup_table);
RcppExport SEXP _utils4ni_count_elegible(SEXP imageSEXP, SEXP patch_sizeSEXP, SEXP search_sizeSEXP, SEXP strideSEXP, SEXP voxel_lookup_tableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type image(imageSEXP);
    Rcpp::traits::input_parameter< int >::type patch_size(patch_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type search_size(search_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type stride(strideSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    rcpp_result_gen = Rcpp::wrap(count_elegible(image, patch_size, search_size, stride, voxel_lookup_table));
    return rcpp_result_gen;
END_RCPP
}
// constrained_initialization_omp
void constrained_initialization_omp(NumericVector input_image, NumericVector template4D, int patch_size, int search_size, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector kANN, int k, int ncores);
RcppExport SEXP _utils4ni_constrained_initialization_omp(SEXP input_imageSEXP, SEXP template4DSEXP, SEXP patch_sizeSEXP, SEXP search_sizeSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP kANNSEXP, SEXP kSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input_image(input_imageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type template4D(template4DSEXP);
    Rcpp::traits::input_parameter< int >::type patch_size(patch_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type search_size(search_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    constrained_initialization_omp(input_image, template4D, patch_size, search_size, actual_voxels, voxel_lookup_table, kANN, k, ncores);
    return R_NilValue;
END_RCPP
}
// all_patches_similarity_omp
void all_patches_similarity_omp(NumericVector input_image, NumericVector template4D, int k, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector patch_neighbours, IntegerVector kANN, NumericVector similarities, int ncores);
RcppExport SEXP _utils4ni_all_patches_similarity_omp(SEXP input_imageSEXP, SEXP template4DSEXP, SEXP kSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP patch_neighboursSEXP, SEXP kANNSEXP, SEXP similaritiesSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input_image(input_imageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type template4D(template4DSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type patch_neighbours(patch_neighboursSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type similarities(similaritiesSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    all_patches_similarity_omp(input_image, template4D, k, actual_voxels, voxel_lookup_table, patch_neighbours, kANN, similarities, ncores);
    return R_NilValue;
END_RCPP
}
// all_patches_similarity_omp2
void all_patches_similarity_omp2(NumericVector input_image, NumericVector template4D, int k, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector voxel_array_index, IntegerVector patch_neighbours, IntegerVector kANN, NumericVector similarities, int ncores);
RcppExport SEXP _utils4ni_all_patches_similarity_omp2(SEXP input_imageSEXP, SEXP template4DSEXP, SEXP kSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP voxel_array_indexSEXP, SEXP patch_neighboursSEXP, SEXP kANNSEXP, SEXP similaritiesSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input_image(input_imageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type template4D(template4DSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_array_index(voxel_array_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type patch_neighbours(patch_neighboursSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type similarities(similaritiesSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    all_patches_similarity_omp2(input_image, template4D, k, actual_voxels, voxel_lookup_table, voxel_array_index, patch_neighbours, kANN, similarities, ncores);
    return R_NilValue;
END_RCPP
}
// propagation_step_omp
void propagation_step_omp(NumericVector input_image, NumericVector template4D, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector patch_neighbours, IntegerVector kANN, int k, int direction, int patch_size, int stride, NumericVector similarities, int ncores);
RcppExport SEXP _utils4ni_propagation_step_omp(SEXP input_imageSEXP, SEXP template4DSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP patch_neighboursSEXP, SEXP kANNSEXP, SEXP kSEXP, SEXP directionSEXP, SEXP patch_sizeSEXP, SEXP strideSEXP, SEXP similaritiesSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input_image(input_imageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type template4D(template4DSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type patch_neighbours(patch_neighboursSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type direction(directionSEXP);
    Rcpp::traits::input_parameter< int >::type patch_size(patch_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type stride(strideSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type similarities(similaritiesSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    propagation_step_omp(input_image, template4D, actual_voxels, voxel_lookup_table, patch_neighbours, kANN, k, direction, patch_size, stride, similarities, ncores);
    return R_NilValue;
END_RCPP
}
// constrained_random_search_omp
void constrained_random_search_omp(NumericVector input_image, NumericVector template4D, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector kANN, int k, int patch_size, IntegerVector patch_neighbours, int search_size_max, NumericVector similarities, int max_random_neighbours, int ncores);
RcppExport SEXP _utils4ni_constrained_random_search_omp(SEXP input_imageSEXP, SEXP template4DSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP kANNSEXP, SEXP kSEXP, SEXP patch_sizeSEXP, SEXP patch_neighboursSEXP, SEXP search_size_maxSEXP, SEXP similaritiesSEXP, SEXP max_random_neighboursSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type input_image(input_imageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type template4D(template4DSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type patch_size(patch_sizeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type patch_neighbours(patch_neighboursSEXP);
    Rcpp::traits::input_parameter< int >::type search_size_max(search_size_maxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type similarities(similaritiesSEXP);
    Rcpp::traits::input_parameter< int >::type max_random_neighbours(max_random_neighboursSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    constrained_random_search_omp(input_image, template4D, actual_voxels, voxel_lookup_table, kANN, k, patch_size, patch_neighbours, search_size_max, similarities, max_random_neighbours, ncores);
    return R_NilValue;
END_RCPP
}
// label_fusion_omp
void label_fusion_omp(IntegerVector labels4D, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector label_ids, IntegerVector kANN, IntegerVector patch_neighbours, int k, double lambda, double sigma2, NumericVector match, NumericVector new_voting, int ncores);
RcppExport SEXP _utils4ni_label_fusion_omp(SEXP labels4DSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP label_idsSEXP, SEXP kANNSEXP, SEXP patch_neighboursSEXP, SEXP kSEXP, SEXP lambdaSEXP, SEXP sigma2SEXP, SEXP matchSEXP, SEXP new_votingSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labels4D(labels4DSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type label_ids(label_idsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type patch_neighbours(patch_neighboursSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type match(matchSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type new_voting(new_votingSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    label_fusion_omp(labels4D, actual_voxels, voxel_lookup_table, label_ids, kANN, patch_neighbours, k, lambda, sigma2, match, new_voting, ncores);
    return R_NilValue;
END_RCPP
}
// label_fusion2_omp
void label_fusion2_omp(IntegerVector labels4D, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector label_ids, IntegerVector kANN, IntegerVector patch_neighbours, int k, double lambda, double sigma2, NumericVector match, NumericVector new_voting, int ncores);
RcppExport SEXP _utils4ni_label_fusion2_omp(SEXP labels4DSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP label_idsSEXP, SEXP kANNSEXP, SEXP patch_neighboursSEXP, SEXP kSEXP, SEXP lambdaSEXP, SEXP sigma2SEXP, SEXP matchSEXP, SEXP new_votingSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labels4D(labels4DSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type label_ids(label_idsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type patch_neighbours(patch_neighboursSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type match(matchSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type new_voting(new_votingSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    label_fusion2_omp(labels4D, actual_voxels, voxel_lookup_table, label_ids, kANN, patch_neighbours, k, lambda, sigma2, match, new_voting, ncores);
    return R_NilValue;
END_RCPP
}
// label_fusion3_omp
void label_fusion3_omp(IntegerVector labels4D, int actual_voxels, IntegerVector voxel_lookup_table, IntegerVector label_ids, IntegerVector kANN, IntegerVector patch_neighbours, int k, double lambda, double sigma2, NumericVector match, NumericVector new_voting, int ncores);
RcppExport SEXP _utils4ni_label_fusion3_omp(SEXP labels4DSEXP, SEXP actual_voxelsSEXP, SEXP voxel_lookup_tableSEXP, SEXP label_idsSEXP, SEXP kANNSEXP, SEXP patch_neighboursSEXP, SEXP kSEXP, SEXP lambdaSEXP, SEXP sigma2SEXP, SEXP matchSEXP, SEXP new_votingSEXP, SEXP ncoresSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labels4D(labels4DSEXP);
    Rcpp::traits::input_parameter< int >::type actual_voxels(actual_voxelsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type voxel_lookup_table(voxel_lookup_tableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type label_ids(label_idsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kANN(kANNSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type patch_neighbours(patch_neighboursSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type match(matchSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type new_voting(new_votingSEXP);
    Rcpp::traits::input_parameter< int >::type ncores(ncoresSEXP);
    label_fusion3_omp(labels4D, actual_voxels, voxel_lookup_table, label_ids, kANN, patch_neighbours, k, lambda, sigma2, match, new_voting, ncores);
    return R_NilValue;
END_RCPP
}
// regularize
NumericVector regularize(NumericVector image, NumericVector kernel);
RcppExport SEXP _utils4ni_regularize(SEXP imageSEXP, SEXP kernelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type image(imageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kernel(kernelSEXP);
    rcpp_result_gen = Rcpp::wrap(regularize(image, kernel));
    return rcpp_result_gen;
END_RCPP
}
// sum_by_ROI
NumericVector sum_by_ROI(IntegerVector labelled, NumericVector values);
RcppExport SEXP _utils4ni_sum_by_ROI(SEXP labelledSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labelled(labelledSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_by_ROI(labelled, values));
    return rcpp_result_gen;
END_RCPP
}
// max_by_ROI
NumericVector max_by_ROI(IntegerVector labelled, NumericVector values);
RcppExport SEXP _utils4ni_max_by_ROI(SEXP labelledSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labelled(labelledSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(max_by_ROI(labelled, values));
    return rcpp_result_gen;
END_RCPP
}
// min_by_ROI
NumericVector min_by_ROI(IntegerVector labelled, NumericVector values);
RcppExport SEXP _utils4ni_min_by_ROI(SEXP labelledSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labelled(labelledSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(min_by_ROI(labelled, values));
    return rcpp_result_gen;
END_RCPP
}
// count_by_ROI
IntegerVector count_by_ROI(IntegerVector labelled);
RcppExport SEXP _utils4ni_count_by_ROI(SEXP labelledSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labelled(labelledSEXP);
    rcpp_result_gen = Rcpp::wrap(count_by_ROI(labelled));
    return rcpp_result_gen;
END_RCPP
}
// mean_by_ROI
NumericVector mean_by_ROI(IntegerVector labelled, NumericVector values);
RcppExport SEXP _utils4ni_mean_by_ROI(SEXP labelledSEXP, SEXP valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type labelled(labelledSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type values(valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_by_ROI(labelled, values));
    return rcpp_result_gen;
END_RCPP
}
// segmentation
NumericVector segmentation(NumericVector image, NumericVector otsu_estimates);
RcppExport SEXP _utils4ni_segmentation(SEXP imageSEXP, SEXP otsu_estimatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type image(imageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type otsu_estimates(otsu_estimatesSEXP);
    rcpp_result_gen = Rcpp::wrap(segmentation(image, otsu_estimates));
    return rcpp_result_gen;
END_RCPP
}
// generate_random
IntegerVector generate_random(int min, int max, int n);
RcppExport SEXP _utils4ni_generate_random(SEXP minSEXP, SEXP maxSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type min(minSEXP);
    Rcpp::traits::input_parameter< int >::type max(maxSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(generate_random(min, max, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_utils4ni_adaptive_nonlocal_means_denoising", (DL_FUNC) &_utils4ni_adaptive_nonlocal_means_denoising, 5},
    {"_utils4ni_cCreateMask", (DL_FUNC) &_utils4ni_cCreateMask, 4},
    {"_utils4ni_cgradientmodule", (DL_FUNC) &_utils4ni_cgradientmodule, 2},
    {"_utils4ni_connected_components", (DL_FUNC) &_utils4ni_connected_components, 1},
    {"_utils4ni_cTruncado", (DL_FUNC) &_utils4ni_cTruncado, 2},
    {"_utils4ni_cvolfilter2d", (DL_FUNC) &_utils4ni_cvolfilter2d, 3},
    {"_utils4ni_defuzzify", (DL_FUNC) &_utils4ni_defuzzify, 1},
    {"_utils4ni_extend_labels", (DL_FUNC) &_utils4ni_extend_labels, 2},
    {"_utils4ni_map_ids_workhorse", (DL_FUNC) &_utils4ni_map_ids_workhorse, 3},
    {"_utils4ni_map_extra_classes", (DL_FUNC) &_utils4ni_map_extra_classes, 3},
    {"_utils4ni_mask_values", (DL_FUNC) &_utils4ni_mask_values, 5},
    {"_utils4ni_get_neighbours", (DL_FUNC) &_utils4ni_get_neighbours, 2},
    {"_utils4ni_fast_read_nifti", (DL_FUNC) &_utils4ni_fast_read_nifti, 1},
    {"_utils4ni_count_elegible", (DL_FUNC) &_utils4ni_count_elegible, 5},
    {"_utils4ni_constrained_initialization_omp", (DL_FUNC) &_utils4ni_constrained_initialization_omp, 9},
    {"_utils4ni_all_patches_similarity_omp", (DL_FUNC) &_utils4ni_all_patches_similarity_omp, 9},
    {"_utils4ni_all_patches_similarity_omp2", (DL_FUNC) &_utils4ni_all_patches_similarity_omp2, 10},
    {"_utils4ni_propagation_step_omp", (DL_FUNC) &_utils4ni_propagation_step_omp, 12},
    {"_utils4ni_constrained_random_search_omp", (DL_FUNC) &_utils4ni_constrained_random_search_omp, 12},
    {"_utils4ni_label_fusion_omp", (DL_FUNC) &_utils4ni_label_fusion_omp, 12},
    {"_utils4ni_label_fusion2_omp", (DL_FUNC) &_utils4ni_label_fusion2_omp, 12},
    {"_utils4ni_label_fusion3_omp", (DL_FUNC) &_utils4ni_label_fusion3_omp, 12},
    {"_utils4ni_regularize", (DL_FUNC) &_utils4ni_regularize, 2},
    {"_utils4ni_sum_by_ROI", (DL_FUNC) &_utils4ni_sum_by_ROI, 2},
    {"_utils4ni_max_by_ROI", (DL_FUNC) &_utils4ni_max_by_ROI, 2},
    {"_utils4ni_min_by_ROI", (DL_FUNC) &_utils4ni_min_by_ROI, 2},
    {"_utils4ni_count_by_ROI", (DL_FUNC) &_utils4ni_count_by_ROI, 1},
    {"_utils4ni_mean_by_ROI", (DL_FUNC) &_utils4ni_mean_by_ROI, 2},
    {"_utils4ni_segmentation", (DL_FUNC) &_utils4ni_segmentation, 2},
    {"_utils4ni_generate_random", (DL_FUNC) &_utils4ni_generate_random, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_utils4ni(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
