


class ErrorTracker {

    constructor(errorDivID, successMsg, callback=null) {
        this.hasError = false;
        this.message = null;
        this.successMsg = successMsg;
        this.target = errorDivID;
        this.callback = callback;
    }

    setTarget(errorDivID) {
        this.target = errorDivID;
    }

    setCallback(callback) {
        this.callback = callback;
    }

    error(message) {
        this.hasError = true;
        this.message = message;

        const target = document.querySelector(`#${this.target}`);
        target.innerHTML = this.message;
        target.style.color = "red";
        target.style.display = "block";

        if (this.callback !== null) {
            this.callback();
        }
    }

    clear() {
        this.hasError = false;
        this.message = this.successMsg;

        const target = document.querySelector(`#${this.target}`);
        target.innerHTML = this.message;
        target.style.color = "green";
        target.style.display = "block";
    }

}


const tracker = new ErrorTracker(null, "Parsing successful!");