art_make(
    LIBRARY_NAME
        ubcc1pi_Analyzers
    LIB_LIBRARIES
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_BASIC}
        ${ROOT_BASIC_LIB_LIST}
        canvas
        cetlib cetlib_except
    MODULE_LIBRARIES
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_BASIC}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
        ${ROOT_BASIC_LIB_LIST}
        canvas
        cetlib cetlib_except
        ubcc1pi_Helpers
    )

install_headers()
install_source()
install_fhicl()
