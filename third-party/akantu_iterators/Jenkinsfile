pipeline {
  parameters {string(defaultValue: '', description: 'api-token', name: 'API_TOKEN')
              string(defaultValue: '', description: 'buildable phid', name: 'BUILD_TARGET_PHID')
              string(defaultValue: '', description: 'Commit id', name: 'COMMIT_ID')
              string(defaultValue: '', description: 'Diff id', name: 'DIFF_ID')
	      string(defaultValue: 'PHID-PROJ-5eqyu6ooyjktagbhf473', description: 'ID of the project', name: 'PROJECT_ID')
  }

  options {
    disableConcurrentBuilds()
  }

  environment {
    PHABRICATOR_HOST = 'https://c4science.ch/api/'
    PYTHONPATH = sh returnStdout: true, script: 'echo ${WORKSPACE}/test/ci/script/'
  }
  
  agent {
    dockerfile {
      filename 'Dockerfile'
      dir 'test/ci'
      additionalBuildArgs '--tag akantu-iterators'
    }
  }
  
  stages {
    stage('Checkout proper commit') {
      steps {
	checkout scm:  [$class: 'GitSCM',
			branches: [[name: "${COMMIT_ID}" ]]
	], changelog: true
      }
    }
        
    stage('Lint') {
      steps {
	sh """
           arc lint --output json --rev HEAD^ | jq . -srM | tee lint.json
           ./test/ci/scripts/hbm send-arc-lint -f lint.json
           """
      }
    }
    
    stage('Configure') {
      steps {
        sh """#!/bin/bash
           set -o pipefail
           mkdir -p build
           cd build
           cmake -DAKANTU_ITERATORS_TESTS:BOOL=ON .. | tee configure.txt
           """
      }
      post {
	failure {
	  uploadArtifact('configure.txt', 'Configure')
	  deleteDir()
	}
      }
    }
    
    stage('Compile tests') {
      steps {
        sh '''#!/bin/bash
           set -o pipefail
	   
           make -C build | tee compilation_test.txt
           '''
      }
      post {
	failure {
	  uploadArtifact('compilation_test.txt', 'Compilation_Tests')
	}
      }
    }

    stage('Tests') {
      steps {
        sh '''
          cd build/
          ctest -T test --no-compress-output || true
          tag=$(head -n 1 < Testing/TAG)
          if [ -e Testing/${tag}/Test.xml ]; then
            cp Testing/${tag}/Test.xml ../CTestResults.xml
          fi
        '''
      }
    }
  }
  post {
    always {
      createArtifact("./CTestResults.xml")

      step([$class: 'XUnitBuilder',
	    thresholds: [
          [$class: 'SkippedThreshold', failureThreshold: '0'],
          [$class: 'FailedThreshold', failureThreshold: '0']],
	    tools: [
	  [$class: 'CTestType', pattern: 'CTestResults.xml', skipNoTestFiles: true]
	]])
    }

    success {
      passed()
    }

    failure {
      failed()
    }
  }
}

def failed() {
  sh "./test/ci/scripts/hbm failed"
}

def passed() {
  sh "./test/ci/scripts/hbm passed"
}

def createArtifact(filename) {
  sh "./test/ci/scripts/hbm send-uri -k 'Jenkins URI' -u ${BUILD_URL} -l 'View Jenkins result'"
  sh "./test/ci/scripts/hbm send-ctest-results -f ${filename}"
}

def uploadArtifact(artifact, name) {
  sh "./test/ci/scripts/hbm upload-file -f ${artifact} -n \"${name}\" -v ${PROJECT_ID}"
}
