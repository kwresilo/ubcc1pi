install_headers()
install_fhicl()
install_source()

file(GLOB xml_file *.xml)
install_fw(LIST ${xml_file})
