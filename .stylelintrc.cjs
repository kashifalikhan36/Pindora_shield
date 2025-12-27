module.exports = {
  extends: ["stylelint-config-tailwindcss"],
  rules: {
    // Disable 'at-rule-no-unknown' because Tailwind introduces several at-rules
    "at-rule-no-unknown": null
  }
};
