pipeline {
    agent {
	dockerfile { additionalBuildArgs '--tag akantu-environment'}
    }

    environment {
	BLA_VENDOR = 'OpenBLAS'
	OMPI_MCA_plm = 'isolated'
	OMPI_MCA_btl = 'self'
    }

    stages {
	stage('Configure') {
	    steps {
		sh 'env'
		sh 'mkdir -p build'
		sh 'cd build; cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE -DAKANTU_IMPLICIT:BOOL=TRUE -DAKANTU_PARALLEL:BOOL=TRUE -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE -DAKANTU_TESTS:BOOL=TRUE ..'
	    }
	}

	stage('Compile') {
	    steps {
		sh 'cd build/src; make || true'
	    }
	}

	stage('Compile tests') {
	    steps {
		sh 'cd build/test; make || true'
	    }
	}
    }
    
    post {
	always {
	    deleteDir()
	}
    }
}
