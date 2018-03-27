pipeline {
    agent {
	dockerfile { additionalBuildArgs '--tag akantu-environment'}
    }

    environment {
	BLA_VENDOR = 'OpenBLAS'
    }

    stages {
	stage('Configure') {
	    steps {
		sh 'mkdir -p build'
		sh 'cd build; cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE -DAKANTU_IMPLICIT:BOOL=TRUE -DAKANTU_PARALLEL:BOOL=TRUE -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE ..'
	    }
	}

	stage('Compile') {
	    steps {
		sh 'cd build; make'
	    }
	}

	stage('Compile tests') {
	    steps {
		sh 'cd build; cmake -DAKANTU_TESTS:BOOL=TRUE ..'
		sh 'cd build; make'
	    }
	}
    }
}
