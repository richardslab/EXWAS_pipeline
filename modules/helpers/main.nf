process checkEnv {
  output:
  stdout

  script:
  """
  python ${baseDir}/modules/helpers/runtime_checks.py
  """
}