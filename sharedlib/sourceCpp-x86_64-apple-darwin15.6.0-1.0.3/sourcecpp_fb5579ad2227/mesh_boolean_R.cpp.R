`.sourceCpp_1_DLLInfo` <- dyn.load('/Users/kyemyungpark/Dropbox/Codes/project_tcell_activation/analysis/density_moving_window/sharedlib/sourceCpp-x86_64-apple-darwin15.6.0-1.0.3/sourcecpp_fb5579ad2227/sourceCpp_2.so')

mesh_boolean <- Rcpp:::sourceCppFunction(function(VAi, FAi, VBi, FBi, op_type) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_mesh_boolean')

rm(`.sourceCpp_1_DLLInfo`)
