`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/kyemyungpark/Dropbox/Codes/project_tcell_activation/analysis/density_moving_window/path_to_current_dir/sharedlib/sourceCpp-x86_64-apple-darwin15.6.0-1.0.3/sourcecpp_164231a1984e2/sourceCpp_4.so')

mesh_boolean <- Rcpp:::sourceCppFunction(function(VAi, FAi, VBi, FBi, op_type) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_mesh_boolean')

rm(`.sourceCpp_1_DLLInfo`)
