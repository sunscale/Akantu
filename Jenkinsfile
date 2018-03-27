pipeline {
  agent {
    dockerfile true
  }
  stages {
    stage('Configure') {
      steps {
        sh 'mkdir build && cd build'
        sh 'cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE -DAKANTU_IMPLICIT:BOOL=TRUE -DAKANTU_PARALLEL:BOOL=TRUE -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE ..'
      }
    }

    stage('Compile') {
      steps {
        sh 'cd build'
        sh 'make'
      }
    }

    stage('Compile tests') {
      steps {
        sh 'cd build'
        sh 'cmake -DAKANTU_TESTS:BOOL=TRUE ..'
	sh 'make'
      }
    }

  }
}
