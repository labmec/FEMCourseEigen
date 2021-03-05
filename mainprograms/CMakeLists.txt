find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(TestEigen mainEigenTest.cpp ${SOURCES2021})
target_link_libraries (TestEigen Eigen3::Eigen ${Boost_LIBRARIES})

add_executable(Fem2021 main2021.cpp ${SOURCES2021})
target_link_libraries (Fem2021 Eigen3::Eigen)

add_executable(TestIntegrate TestIntegrate.cpp ${SOURCES2021})
target_link_libraries (TestIntegrate Eigen3::Eigen)

add_executable(TestOneDProblem TestOneDProblem.cpp ${SOURCES2021})
target_link_libraries (TestOneDProblem Eigen3::Eigen)

add_executable(TestTwoDProblem TestTwoDProblem.cpp ${SOURCES2021})
target_link_libraries (TestTwoDProblem Eigen3::Eigen)