import java.util.regex.Pattern

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
	stage('CheckOut') {
	    steps {
		checkout scm
	    }
	}
	stage('BuildAndTest') {
    matrix {
	axes {
	    axis {
		name 'NSPACE'
		values 'real', 'complex'
	    }
	    axis {
		name 'PRECISION'
		values 'full', 'mixed'
	    }
	}
    stages {
        stage('Build') {
            steps {
		echo 'building...'
		dir ('./build')
		{
		    sh "../tests/test_automation/jenkins_build_cpu.sh ${NSPACE} ${PRECISION}"
		}
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
		dir('./build')
		{
		    sh '../test/test_automation/jenkins_test.sh ${NSPACE} ${PRECISION}'
		}
            }
	}
	}
    }
}
    }
}
