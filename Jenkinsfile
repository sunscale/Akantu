pipeline {
  parameters {string(defaultValue: '', description: 'api-token', name: 'API_TOKEN')
              string(defaultValue: '', description: 'buildable phid', name: 'BUILD_TARGET_PHID')
              string(defaultValue: '', description: 'Commit id', name: 'COMMIT_ID')
              string(defaultValue: '', description: 'Diff id', name: 'DIFF_ID')
  }

  options {
    disableConcurrentBuilds()
  }

  environment {
    PHABRICATOR_HOST = 'https://c4science.ch/api/'
    PYTHONPATH = sh returnStdout: true, script: 'echo ${WORKSPACE}/test/ci/script/'
    BLA_VENDOR = 'OpenBLAS'
    OMPI_MCA_plm = 'isolated'
    OMPI_MCA_btl = 'tcp,self'
  }
  
  agent {
    dockerfile {
      filename 'Dockerfile'
      dir 'test/ci'
      additionalBuildArgs '--tag akantu-environment'
    }
  }
  
  stages {
    stage('Lint') {
      steps {
	script {
	  try {
	    sh """
               arc lint --output json --rev HEAD^${}^ | jq . -srM | tee lint.json
               if [ -n "$(grep '\[\]' lint.json)" ]; then
                  exit 1
               fi
               """
	  }
	  catch (exc) {
      	    sh """
               ./test/ci/scripts/hbm send-arc-lint -f lint.json
               rm lint.json
               """
	    currentBuild.result = 'UNSTABLE'
	  }
	}
      }
    }
    stage('Configure') {
      steps {
        sh """
        mkdir -p build
        cd build
        cmake -DAKANTU_COHESIVE_ELEMENT:BOOL=TRUE \
              -DAKANTU_IMPLICIT:BOOL=TRUE \
              -DAKANTU_PARALLEL:BOOL=TRUE \
              -DAKANTU_PYTHON_INTERFACE:BOOL=TRUE \
              -DAKANTU_TESTS:BOOL=TRUE ..
        """
      }
      post {
	failure {
	  deleteDir()
	}
      }
    }
    stage('Compile') {
      steps {
	sh 'make -C build/src || true'
      }
    }

    stage ('Warnings gcc') {
      steps {
        warnings(consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']])
      }
    }

    stage('Compile python') {
      steps {
        sh 'make -C build/python || true'
      }
    }

    stage('Compile tests') {
      steps {
        sh 'make -C build/test || true'
      }
    }

    stage('Tests') {
      steps {
        sh '''
          rm -rf build/gtest_reports
          cd build/
          #source ./akantu_environement.sh
        
          ctest -T test --no-compress-output || true
        '''
      }
      post {
	always {
	  script {
	    def TAG = sh returnStdout: true, script: 'head -n 1 < build/Testing/TAG'
	    def TAG_ = TAG.trim()

	    if (fileExists("build/Testing/${TAG}/Test.xml")) {
	      sh "cp build/Testing/${TAG}/Test.xml CTestResults.xml"
	    }
	  }
	}
      }
    }
  }

  post {
    always {
      step([$class: 'XUnitBuilder',
	    thresholds: [
          [$class: 'SkippedThreshold', failureThreshold: '0'],
          [$class: 'FailedThreshold', failureThreshold: '0']],
	    tools: [
	  [$class: 'CTestType', pattern: 'CTestResults.xml', skipNoTestFiles: true]
	]])
      step([$class: 'XUnitBuilder',
            thresholds: [
          [$class: 'SkippedThreshold', failureThreshold: '100'],
          [$class: 'FailedThreshold', failureThreshold: '0']],
            tools: [
	  [$class: 'GoogleTestType', pattern: 'build/gtest_reports/**', skipNoTestFiles: true]
	]])

      createArtifact("./CTestResults.xml")     
    }

    success {
      passed()
    }

    failure {
      emailext(
        body: '''${SCRIPT, template="groovy-html.template"}''',
	mimeType: 'text/html',
        subject: "[Jenkins] ${currentBuild.fullDisplayName} Failed",
	recipientProviders: [[$class: 'CulpritsRecipientProvider']],
	to: 'akantu-admins@akantu.ch',
	replyTo: 'akantu-admins@akantu.ch',
	attachLog: true,
        compressLog: false)
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

def createArtifact(artefact) {
  sh "./test/ci/scripts/hbm send-uri -k 'Jenkins URI' -u ${BUILD_URL} -l 'View Jenkins result'"
  sh "./test/ci/scripts/hbm send-ctest-results -f ${artefact}"
}
