pipeline {
    agent {
	any {
	    node {
		label 'oxygen CI'
		customWorkspace '/dev/shm/jenkins'
	    }
	}
    }
    environment {
	LD_LIBRARY_PATH="""${sh(
returnStdout: true,
script: 'echo "/usr/local/lib:${LD_LIBRARY_PATH}"'
)}"""
    }
    options {
	buildDiscarder(logRotator(numToKeepStr: '10'))
    }
    stages {
        stage('Build') {
            steps {
		echo 'building...'
		checkout scm
		dir './build'
		{
		    sh 'cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DENABLE_SOA=0 -DBUILD_AFQMC=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 .. 2>&1 | tee cmake.out'
		    sh 'make -j12'
		}
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
		dir './build'
		{
		    sh 'ctest -L unit --output-on-failure --timeout 120'
		}
            }
        }
    }
}

