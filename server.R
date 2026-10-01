function(input, output, session) {
  result <- mod_data_server("data")
  mod_qc_server("qc", result)
  mod_norm_server("norm", result)
  mod_de_server("de", result)
}
