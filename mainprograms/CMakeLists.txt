find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(TestEigen mainEigenTest.cpp ${FemSCSOURCE})
target_link_libraries (TestEigen Eigen3::Eigen)

add_executable(Fem2021 main2021.cpp ${FemSOURCE2021} ${SOURCES2021})
target_link_libraries (Fem2021 Eigen3::Eigen)

# NEED TO CHECK THESE TARGETS
# add_executable(TestIntegrate TestIntegrate.cpp ${FemSCSOURCE})

# add_executable(TestOneDProblem TestOneDProblem.cpp ${FemSCSOURCE})

# add_executable(TestTwoDProblem TestTwoDProblem.cpp ${FemSCSOURCE})

# add_executable(TesMatrixXd TesMatrixXd.cpp ${FemSCSOURCE})
# target_link_libraries (TesMatrixXd Eigen3::Eigen)

