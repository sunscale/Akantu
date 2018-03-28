pipeline {
  agent {
    dockerfile {
      additionalBuildArgs '--tag akantu-environment'
    }
    
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
    stage('') {
      steps {
        emailext(subject: 'Build failed for job ${currentBuild.fullDisplayName}', to: 'akantu-admins@akantu.ch', attachLog: true, compressLog: true, body: 'Something is wrong with ${env.BUILD_URL}, failed at step')
      }
    }
  }
  environment {
    BLA_VENDOR = 'OpenBLAS'
    OMPI_MCA_plm = 'isolated'
    OMPI_MCA_btl = 'self'
  }
  post {
    always {
      deleteDir()
      
    }
    
  }
}