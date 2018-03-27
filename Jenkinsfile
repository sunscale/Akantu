pipeline {
  agent {
    docker {
      image 'richart/akantu-public-debian'
    }
  }
  stages {
    stage('Configure') {
      steps {
        sh 'mkdir build && cd build'
        sh 'cmake ..'
      }
    }
  }
}
