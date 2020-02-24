pipeline {
    agent {
	node {
	    label 'master'
	    customWorkspace '/dev/shm/jenkins'
	}
    }
    environment {
	LD_LIBRARY_PATH = """${sh(
returnStdout: true,
script: 'echo "/data/epd/spack/opt/spack/linux-rhel7-x86_64/gcc-7.4.0/openmpi-3.1.3-lzlzqpa4gfnarehk2knpe4fbm3xujstc/lib:/data/epd/spack/opt/spack/linux-rhel7-x86_64/gcc-7.3.0/hdf5-1.10.4-gg33hqk5cbr4kmza3hbal2wa7bv5hrh2/lib:/usr/local/lib:/data/epd/spack/opt/spack/linux-rhel7-x86_64/gcc-7.3.0/boost-1.69.0-ffo2xwbslhmjfg4v6cthxs5ypz6pq5em/lib:${LD_LIBRARY_PATH}"'
)}"""
	PATH = """${sh(
returnStdout: true,
script: 'echo "/data/epd/spack/opt/spack/linux-rhel7-x86_64/gcc-7.4.0/openmpi-3.1.3-lzlzqpa4gfnarehk2knpe4fbm3xujstc/bin:/data/epd/spack/opt/spack/linux-rhel7-x86_64/gcc-7.3.0/hdf5-1.10.4-gg33hqk5cbr4kmza3hbal2wa7bv5hrh2/bin:${PATH}"'
)}"""
	CMAKE_PREFIX_PATH = """${sh(
returnStdout: true,
script: 'echo "/data/epd/spack/opt/spack/linux-rhel7-x86_64/gcc-7.3.0/boost-1.69.0-ffo2xwbslhmjfg4v6cthxs5ypz6pq5em/:${CMAKE_PREFIX_PATH}"'
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
		dir ('./build')
		{
		    sh 'cmake -DQMC_COMPLEX=0 -DQMC_MIXED_PRECISION=0 -DENABLE_SOA=0 -DBUILD_AFQMC=1 -DCMAKE_C_COMPILER="mpicc" -DCMAKE_CXX_COMPILER="mpicxx" -DQMC_NO_SLOW_CUSTOM_TESTING_COMMANDS=1 .. 2>&1 | tee cmake.out'
		    sh 'make -j12'
		}
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
		dir('./build')
		{
		    sh 'ctest -L unit --output-on-failure --timeout 120'
		}
            }
	}
    }
}

