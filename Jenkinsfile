Jenkinsfile (Declarative Pipeline)

pipeline {
    agent any
    environment {
	LD_LIBRARY_PATH="""${sh(
returnStdout: true,
script: 'echo "/usr/local/lib:${LD_LIBRARY_PATH}"'
)}"""
    }
    stages {
        stage('Build') {
            steps {
		echo 'building...'
		checkout scm
		sh 'cd qmcpack'
		sh 'cd build'
		sh 'cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DENABLE_SOA=0 -DBUILD_AFQMC=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 .. 2>&1 | tee cmake.out'
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
		sh 'ctest -L unit --output-on-failure --timeout 120'
            }
        }
    }
}

